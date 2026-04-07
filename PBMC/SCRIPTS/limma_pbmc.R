################################################################################
# STAGE 10: CELL TYPE COMPOSITION + DEG (LIMMA / VOOM)
################################################################################

stage10_composition_and_deg <- function(
    obj,
    out_dir = file.path(OUT_DIR, "stage10_limma_deg"),
    celltype_col = NULL,
    assay = "RNA",
    min_cells_per_pb = 20,
    min_samples_per_group = 2,
    use_cd45_only = FALSE
) {
  section_header("STAGE 10: CELL TYPE COMPOSITION + DEG (LIMMA)")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --------------------------------------------------------------------------
  # Helpers
  # --------------------------------------------------------------------------
  .get_counts_matrix <- function(seu, assay = "RNA") {
    seu <- join_if_multi(seu, assay = assay, stage_name = "JOIN_BEFORE_STAGE10")
    
    # Seurat v5-safe counts extraction
    counts_layers <- grep("^counts", SeuratObject::Layers(seu[[assay]]), value = TRUE)
    
    if (length(counts_layers) > 0) {
      mat <- SeuratObject::LayerData(seu[[assay]], layer = counts_layers[1])
    } else {
      mat <- SeuratObject::GetAssayData(seu, assay = assay, slot = "counts")
    }
    
    if (!inherits(mat, "dgCMatrix")) {
      mat <- as(mat, "dgCMatrix")
    }
    return(mat)
  }
  
  .safe_topTable <- function(fit, coef_name) {
    tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
    tt$gene <- rownames(tt)
    tt
  }
  
  .write_skip <- function(path, reason) {
    data.table::fwrite(data.frame(reason = reason), path)
  }
  
  # --------------------------------------------------------------------------
  # Metadata setup
  # --------------------------------------------------------------------------
  md <- obj@meta.data |>
    tibble::rownames_to_column("cell") |>
    as.data.frame(stringsAsFactors = FALSE)
  
  if (is.null(celltype_col)) {
    candidate_cols <- c("Monaco_main_pruned", "Monaco_main", "cluster_annotation", "cluster")
    celltype_col <- candidate_cols[candidate_cols %in% colnames(md)][1]
  }
  
  if (length(celltype_col) == 0 || is.na(celltype_col) || !(celltype_col %in% colnames(md))) {
    stop("No valid celltype column found. Tried: Monaco_main_pruned, Monaco_main, cluster_annotation, cluster")
  }
  
  log_status("STAGE10", paste("Using cell type column:", celltype_col), "INFO")
  
  md <- md |>
    dplyr::mutate(
      subject_id = dplyr::case_when(
        !is.na(Subject_raw) & Subject_raw != "" ~ as.character(Subject_raw),
        !is.na(Subject_fixed) & Subject_fixed != "" ~ as.character(Subject_fixed),
        TRUE ~ NA_character_
      ),
      group = dplyr::case_when(
        Treatment == "HC"        ~ "HC",
        Treatment == "Untreated" ~ "Mild_Untreated",
        Treatment == "Treated"   ~ "Mild_Treated",
        Treatment == "Severe"    ~ "Severe",
        TRUE ~ NA_character_
      ),
      celltype = .data[[celltype_col]]
    ) |>
    dplyr::mutate(
      celltype = ifelse(is.na(celltype) | celltype == "", "Unassigned", as.character(celltype)),
      # safe sample IDs that won't get rewritten
      sample_unit = dplyr::case_when(
        group == "HC"             ~ paste0(subject_id, ".HC"),
        group == "Severe"         ~ paste0(subject_id, ".SEV"),
        group == "Mild_Untreated" ~ paste0(subject_id, ".UNT"),
        group == "Mild_Treated"   ~ paste0(subject_id, ".TRT"),
        TRUE ~ NA_character_
      )
    )
  
  if (use_cd45_only && "CD45pos" %in% colnames(md)) {
    md <- md |>
      dplyr::filter(CD45pos %in% TRUE)
    log_status("STAGE10", "Restricted to CD45+ cells only", "INFO")
  }
  
  md <- md |>
    dplyr::filter(!is.na(cell), !is.na(group), !is.na(subject_id), !is.na(sample_unit), !is.na(celltype))
  
  data.table::fwrite(md, file.path(out_dir, "stage10_cell_metadata_used.csv"))
  
  sample_meta <- md |>
    dplyr::distinct(sample_unit, subject_id, group) |>
    dplyr::arrange(group, subject_id)
  
  data.table::fwrite(sample_meta, file.path(out_dir, "stage10_sample_metadata.csv"))
  
  # --------------------------------------------------------------------------
  # PART A: Cell type composition
  # --------------------------------------------------------------------------
  log_status("COMPOSITION", "Running cell type composition analysis", "PROGRESS")
  
  comp_counts <- md |>
    dplyr::count(sample_unit, subject_id, group, celltype, name = "n_cells")
  
  total_cells <- comp_counts |>
    dplyr::group_by(sample_unit) |>
    dplyr::summarise(total_cells = sum(n_cells), .groups = "drop")
  
  comp_counts <- comp_counts |>
    dplyr::left_join(total_cells, by = "sample_unit") |>
    dplyr::mutate(
      prop = n_cells / total_cells,
      logit_prop = qlogis((n_cells + 0.5) / (total_cells + 1))
    )
  
  data.table::fwrite(comp_counts, file.path(out_dir, "composition_sample_by_celltype_counts.csv"))
  
  comp_mat <- comp_counts |>
    dplyr::select(celltype, sample_unit, logit_prop) |>
    tidyr::pivot_wider(names_from = sample_unit, values_from = logit_prop, values_fill = 0)
  
  comp_mat2 <- as.data.frame(comp_mat)
  rownames(comp_mat2) <- comp_mat2$celltype
  comp_mat2$celltype <- NULL
  comp_mat2 <- as.matrix(comp_mat2)
  
  # 3-group composition
  comp_meta_3g <- sample_meta |>
    dplyr::filter(group %in% c("HC", "Mild_Untreated", "Severe")) |>
    dplyr::arrange(match(group, c("Mild_Untreated", "HC", "Severe")), subject_id)
  
  if (nrow(comp_meta_3g) > 0) {
    mat_3g <- comp_mat2[, comp_meta_3g$sample_unit, drop = FALSE]
    comp_meta_3g$group <- factor(comp_meta_3g$group, levels = c("Mild_Untreated", "HC", "Severe"))
    
    design_3g <- model.matrix(~0 + group, data = comp_meta_3g)
    colnames(design_3g) <- make.names(colnames(design_3g))
    
    fit_3g <- limma::lmFit(mat_3g, design_3g)
    
    cont_3g <- limma::makeContrasts(
      HC_vs_MildUntreated     = groupHC - groupMild_Untreated,
      Severe_vs_MildUntreated = groupSevere - groupMild_Untreated,
      Severe_vs_HC            = groupSevere - groupHC,
      levels = design_3g
    )
    
    fit_3g <- limma::contrasts.fit(fit_3g, cont_3g)
    fit_3g <- limma::eBayes(fit_3g)
    
    for (coef_name in colnames(cont_3g)) {
      tt <- limma::topTable(fit_3g, coef = coef_name, number = Inf, sort.by = "P")
      tt$celltype <- rownames(tt)
      data.table::fwrite(tt, file.path(out_dir, paste0("composition_", coef_name, ".csv")))
    }
  }
  
  # paired composition
  comp_meta_pair <- sample_meta |>
    dplyr::filter(group %in% c("Mild_Untreated", "Mild_Treated")) |>
    dplyr::arrange(subject_id, group)
  
  if (nrow(comp_meta_pair) > 0) {
    pair_subjects <- comp_meta_pair |>
      dplyr::count(subject_id, group) |>
      tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0) |>
      dplyr::filter(Mild_Untreated > 0, Mild_Treated > 0) |>
      dplyr::pull(subject_id)
    
    comp_meta_pair <- comp_meta_pair |>
      dplyr::filter(subject_id %in% pair_subjects)
    
    if (nrow(comp_meta_pair) > 0) {
      mat_pair <- comp_mat2[, comp_meta_pair$sample_unit, drop = FALSE]
      comp_meta_pair$group <- factor(comp_meta_pair$group, levels = c("Mild_Untreated", "Mild_Treated"))
      comp_meta_pair$subject_id <- factor(comp_meta_pair$subject_id)
      
      design_pair <- model.matrix(~ subject_id + group, data = comp_meta_pair)
      fit_pair <- limma::lmFit(mat_pair, design_pair)
      fit_pair <- limma::eBayes(fit_pair)
      
      coef_name <- grep("^group", colnames(design_pair), value = TRUE)[1]
      tt_pair <- limma::topTable(fit_pair, coef = coef_name, number = Inf, sort.by = "P")
      tt_pair$celltype <- rownames(tt_pair)
      data.table::fwrite(tt_pair, file.path(out_dir, "composition_MildTreated_vs_MildUntreated_paired.csv"))
    }
  }
  
  # --------------------------------------------------------------------------
  # PART B: Manual pseudobulk
  # --------------------------------------------------------------------------
  log_status("DEG", "Building pseudobulk counts", "PROGRESS")
  
  counts_mat <- .get_counts_matrix(obj, assay = assay)
  
  common_cells <- intersect(colnames(counts_mat), md$cell)
  counts_mat <- counts_mat[, common_cells, drop = FALSE]
  md <- md |>
    dplyr::filter(cell %in% common_cells)
  
  md <- md |>
    dplyr::mutate(pb_id = paste(celltype, sample_unit, sep = "|||"))
  
  pb_groups <- split(md$cell, md$pb_id)
  
  if (length(pb_groups) == 0) {
    stop("No pseudobulk groups were created.")
  }
  
  pb_counts_list <- lapply(pb_groups, function(cells_use) {
    cells_use <- intersect(cells_use, colnames(counts_mat))
    if (length(cells_use) == 0) {
      return(rep(0, nrow(counts_mat)))
    }
    Matrix::rowSums(counts_mat[, cells_use, drop = FALSE])
  })
  
  agg <- do.call(cbind, pb_counts_list)
  if (is.null(dim(agg))) {
    agg <- matrix(agg, ncol = 1)
  }
  rownames(agg) <- rownames(counts_mat)
  colnames(agg) <- names(pb_groups)
  agg <- as.matrix(agg)
  
  pb_meta <- data.frame(pb_id = colnames(agg), stringsAsFactors = FALSE) |>
    tidyr::separate(pb_id, into = c("celltype", "sample_unit"), sep = "\\|\\|\\|", extra = "merge", remove = FALSE) |>
    dplyr::left_join(sample_meta, by = "sample_unit")
  
  pb_cell_counts <- md |>
    dplyr::count(celltype, sample_unit, name = "n_cells_pb")
  
  pb_meta <- pb_meta |>
    dplyr::left_join(pb_cell_counts, by = c("celltype", "sample_unit")) |>
    dplyr::mutate(n_cells_pb = ifelse(is.na(n_cells_pb), 0L, n_cells_pb)) |>
    dplyr::arrange(celltype, group, subject_id)
  
  data.table::fwrite(pb_meta, file.path(out_dir, "pseudobulk_metadata.csv"))
  
  # --------------------------------------------------------------------------
  # PART C: DEG
  # --------------------------------------------------------------------------
  log_status("DEG", "Running cell-type-specific pseudobulk DEG", "PROGRESS")
  
  run_deg_one_celltype <- function(ct, counts_mat_pb, pb_meta_df, out_dir_ct) {
    dir.create(out_dir_ct, recursive = TRUE, showWarnings = FALSE)
    
    ct_meta <- pb_meta_df |>
      dplyr::filter(celltype == ct, n_cells_pb >= min_cells_per_pb)
    
    if (nrow(ct_meta) == 0) {
      .write_skip(file.path(out_dir_ct, paste0("SKIPPED_", make.names(ct), ".csv")),
                  "No pseudobulk samples passed min_cells_per_pb")
      return(invisible(NULL))
    }
    
    ct_cols <- ct_meta$pb_id[ct_meta$pb_id %in% colnames(counts_mat_pb)]
    if (length(ct_cols) == 0) {
      .write_skip(file.path(out_dir_ct, paste0("SKIPPED_", make.names(ct), ".csv")),
                  "No pseudobulk columns matched this cell type")
      return(invisible(NULL))
    }
    
    ct_counts <- counts_mat_pb[, ct_cols, drop = FALSE]
    
    if (ncol(ct_counts) == 0 || length(ct_counts) == 0) {
      .write_skip(file.path(out_dir_ct, paste0("SKIPPED_", make.names(ct), ".csv")),
                  "Empty pseudobulk count matrix")
      return(invisible(NULL))
    }
    
    # ---------------------- 3-group unpaired ----------------------
    ct_meta_3g <- ct_meta |>
      dplyr::filter(group %in% c("HC", "Mild_Untreated", "Severe"))
    
    grp_n <- table(ct_meta_3g$group)
    
    if (all(c("HC", "Mild_Untreated", "Severe") %in% names(grp_n)) &&
        all(grp_n[c("HC", "Mild_Untreated", "Severe")] >= min_samples_per_group)) {
      
      ct_meta_3g <- ct_meta_3g |>
        dplyr::arrange(match(group, c("Mild_Untreated", "HC", "Severe")), subject_id)
      
      y_counts <- ct_counts[, ct_meta_3g$pb_id, drop = FALSE]
      
      if (ncol(y_counts) > 0 && sum(y_counts) > 0) {
        ct_meta_3g$group <- factor(ct_meta_3g$group, levels = c("Mild_Untreated", "HC", "Severe"))
        design <- model.matrix(~0 + group, data = ct_meta_3g)
        colnames(design) <- make.names(colnames(design))
        
        dge <- edgeR::DGEList(counts = y_counts)
        keep_gene <- edgeR::filterByExpr(dge, design = design)
        dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
        
        if (nrow(dge) > 0) {
          dge <- edgeR::calcNormFactors(dge)
          v <- limma::voom(dge, design, plot = FALSE)
          fit <- limma::lmFit(v, design)
          
          cont <- limma::makeContrasts(
            HC_vs_MildUntreated     = groupHC - groupMild_Untreated,
            Severe_vs_MildUntreated = groupSevere - groupMild_Untreated,
            Severe_vs_HC            = groupSevere - groupHC,
            levels = design
          )
          
          fit2 <- limma::contrasts.fit(fit, cont)
          fit2 <- limma::eBayes(fit2)
          
          for (coef_name in colnames(cont)) {
            tt <- .safe_topTable(fit2, coef_name)
            data.table::fwrite(
              tt,
              file.path(out_dir_ct, paste0("DEG_", make.names(ct), "_", coef_name, ".csv"))
            )
          }
        } else {
          .write_skip(file.path(out_dir_ct, paste0("SKIPPED_", make.names(ct), "_3group.csv")),
                      "No genes passed filterByExpr for 3-group comparison")
        }
      } else {
        .write_skip(file.path(out_dir_ct, paste0("SKIPPED_", make.names(ct), "_3group.csv")),
                    "3-group pseudobulk matrix had zero columns or zero total counts")
      }
    } else {
      .write_skip(file.path(out_dir_ct, paste0("SKIPPED_", make.names(ct), "_3group.csv")),
                  "Insufficient samples for HC vs Mild_Untreated vs Severe")
    }
    
    # ---------------------- paired mild untreated vs treated ----------------------
    ct_meta_pair <- ct_meta |>
      dplyr::filter(group %in% c("Mild_Untreated", "Mild_Treated"))
    
    if (nrow(ct_meta_pair) > 0) {
      pair_subjects <- ct_meta_pair |>
        dplyr::count(subject_id, group) |>
        tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0) |>
        dplyr::filter(Mild_Untreated > 0, Mild_Treated > 0) |>
        dplyr::pull(subject_id)
      
      ct_meta_pair <- ct_meta_pair |>
        dplyr::filter(subject_id %in% pair_subjects) |>
        dplyr::arrange(subject_id, group)
      
      if (length(unique(ct_meta_pair$subject_id)) >= 2) {
        y_counts2 <- ct_counts[, ct_meta_pair$pb_id, drop = FALSE]
        
        if (ncol(y_counts2) > 0 && sum(y_counts2) > 0) {
          ct_meta_pair$group <- factor(ct_meta_pair$group, levels = c("Mild_Untreated", "Mild_Treated"))
          ct_meta_pair$subject_id <- factor(ct_meta_pair$subject_id)
          
          design <- model.matrix(~ subject_id + group, data = ct_meta_pair)
          
          dge <- edgeR::DGEList(counts = y_counts2)
          keep_gene <- edgeR::filterByExpr(dge, design = design)
          dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
          
          if (nrow(dge) > 0) {
            dge <- edgeR::calcNormFactors(dge)
            v <- limma::voom(dge, design, plot = FALSE)
            fit <- limma::lmFit(v, design)
            fit <- limma::eBayes(fit)
            
            coef_name <- grep("^group", colnames(design), value = TRUE)[1]
            
            if (!is.na(coef_name)) {
              tt <- .safe_topTable(fit, coef_name)
              data.table::fwrite(
                tt,
                file.path(out_dir_ct, paste0("DEG_", make.names(ct), "_MildTreated_vs_MildUntreated_paired.csv"))
              )
            } else {
              .write_skip(file.path(out_dir_ct, paste0("SKIPPED_", make.names(ct), "_paired.csv")),
                          "No treatment coefficient found in paired design")
            }
          } else {
            .write_skip(file.path(out_dir_ct, paste0("SKIPPED_", make.names(ct), "_paired.csv")),
                        "No genes passed filterByExpr for paired comparison")
          }
        } else {
          .write_skip(file.path(out_dir_ct, paste0("SKIPPED_", make.names(ct), "_paired.csv")),
                      "Paired pseudobulk matrix had zero columns or zero total counts")
        }
      } else {
        .write_skip(file.path(out_dir_ct, paste0("SKIPPED_", make.names(ct), "_paired.csv")),
                    "Fewer than 2 complete paired mild subjects")
      }
    } else {
      .write_skip(file.path(out_dir_ct, paste0("SKIPPED_", make.names(ct), "_paired.csv")),
                  "No mild untreated/treated samples for paired comparison")
    }
    
    invisible(NULL)
  }
  
  celltypes_to_run <- pb_meta |>
    dplyr::filter(n_cells_pb >= min_cells_per_pb) |>
    dplyr::count(celltype, sort = TRUE) |>
    dplyr::pull(celltype)
  
  data.table::fwrite(
    data.frame(celltype = celltypes_to_run),
    file.path(out_dir, "celltypes_tested_for_deg.csv")
  )
  
  deg_root <- file.path(out_dir, "deg_by_celltype")
  dir.create(deg_root, recursive = TRUE, showWarnings = FALSE)
  
  for (ct in celltypes_to_run) {
    log_status("DEG", paste("Running:", ct), "INFO")
    ct_dir <- file.path(deg_root, make.names(ct))
    run_deg_one_celltype(ct, agg, pb_meta, ct_dir)
  }
  
  log_status("STAGE10", "Stage 10 complete", "SUCCESS")
  
  return(list(
    cell_metadata = md,
    sample_metadata = sample_meta,
    pseudobulk_metadata = pb_meta,
    pseudobulk_counts = agg
  ))
}


stage10_results <- stage10_composition_and_deg(
  obj = obj_filtered,
  out_dir = file.path(OUT_DIR, "stage10_limma_deg"),
  celltype_col = "Monaco_main_pruned",
  assay = "RNA",
  min_cells_per_pb = 20,
  min_samples_per_group = 2,
  use_cd45_only = FALSE
)