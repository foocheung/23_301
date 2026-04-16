###############################################################################
# PSEUDOBULK DEG + GSEA SIDE-BY-SIDE
# Runs:
#   1) pooled unpaired DEG/GSEA
#   2) within-pool unpaired DEG/GSEA
#   3) paired cross-pool DEG/GSEA
#
# GSEA databases:
#   - Hallmark
#   - GO Biological Process
#   - Reactome
#
# Notes for your design:
#   - Mild Untreated is only in Pool1
#   - Mild Treated is only in Pool2
#   - HC is in Pool1 and Pool2
#   - Severe is in Pool1 and Pool2
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(msigdbr)
  library(enrichplot)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(forcats)
})

# =============================================================================
# SETTINGS
# =============================================================================
MIN_CELLS_PER_PB   <- 20
MIN_SAMPLES_GROUP  <- 2

DEG_FDR_CUTOFF     <- 0.05
DEG_LOGFC_CUTOFF   <- 0.25

GSEA_MIN_GS        <- 10
GSEA_MAX_GS        <- 500
GSEA_PADJ_CUTOFF   <- 0.05
GSEA_TOP_N_PLOT    <- 15

# =============================================================================
# HELPERS
# =============================================================================
log_msg <- function(section, msg) {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%H:%M:%S"), section, msg))
}

safe_name <- function(x) {
  x <- gsub("[/\\\\]", "_", x)
  x <- gsub("[^A-Za-z0-9._-]", "_", x)
  x
}

save_plot <- function(p, path_no_ext, width = 10, height = 7) {
  ggplot2::ggsave(paste0(path_no_ext, ".png"), p, width = width, height = height, dpi = 300)
  ggplot2::ggsave(paste0(path_no_ext, ".pdf"), p, width = width, height = height)
}

write_skip <- function(path, reason) {
  data.table::fwrite(data.frame(reason = reason), path)
}

get_counts_matrix_v5 <- function(obj, assay = "RNA") {
  assay_obj <- obj[[assay]]
  
  if ("Assay5" %in% class(assay_obj)) {
    layer_names <- SeuratObject::Layers(assay_obj)
    count_layers <- grep("^counts", layer_names, value = TRUE)
    if (length(count_layers) == 0) {
      stop("No counts layer found in assay: ", assay)
    }
    mat <- SeuratObject::LayerData(assay_obj, layer = count_layers[1])
  } else {
    mat <- SeuratObject::GetAssayData(obj, assay = assay, layer = "counts")
  }
  
  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
  }
  mat
}

symbol_to_entrez_ranked <- function(gene_list_symbol) {
  map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = names(gene_list_symbol),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  )
  
  map <- map %>%
    dplyr::filter(!is.na(SYMBOL), !is.na(ENTREZID)) %>%
    dplyr::distinct(SYMBOL, ENTREZID)
  
  if (nrow(map) == 0) return(NULL)
  
  common_syms <- intersect(names(gene_list_symbol), map$SYMBOL)
  if (length(common_syms) == 0) return(NULL)
  
  map2 <- map %>%
    dplyr::filter(SYMBOL %in% common_syms)
  
  map2$stat <- unname(gene_list_symbol[map2$SYMBOL])
  
  map2 <- map2 %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  map2 <- map2 %>%
    dplyr::group_by(ENTREZID) %>%
    dplyr::slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  out <- map2$stat
  names(out) <- map2$ENTREZID
  
  out <- as.numeric(out)
  names(out) <- map2$ENTREZID
  out <- out[!is.na(out) & !is.na(names(out))]
  out <- sort(out, decreasing = TRUE)
  
  return(out)
}

prep_gsea_plot_df <- function(df, top_n = 15) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df <- as.data.frame(df)
  needed <- c("Description", "NES", "p.adjust")
  if (!all(needed %in% colnames(df))) return(NULL)
  
  df <- df %>%
    dplyr::filter(!is.na(NES), !is.na(p.adjust)) %>%
    dplyr::arrange(p.adjust, dplyr::desc(abs(NES)))
  
  if (nrow(df) == 0) return(NULL)
  
  df <- utils::head(df, min(top_n, nrow(df)))
  df$pathway_short <- stringr::str_trunc(df$Description, width = 60)
  df$pathway_short <- make.unique(df$pathway_short)
  
  ord <- order(df$NES, na.last = TRUE)
  df$pathway_short <- factor(df$pathway_short, levels = df$pathway_short[ord])
  df$direction <- ifelse(df$NES > 0, "Higher_in_first_group", "Higher_in_second_group")
  df
}

plot_gsea_bar <- function(gsea_df, title_txt, out_prefix, top_n = 15) {
  plot_df <- prep_gsea_plot_df(gsea_df, top_n = top_n)
  if (is.null(plot_df) || nrow(plot_df) == 0) return(invisible(NULL))
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = pathway_short, y = NES, fill = NES > 0)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
    ggplot2::labs(title = title_txt, x = NULL, y = "NES") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "none")
  
  save_plot(p, out_prefix, width = 10, height = 7)
  invisible(p)
}

plot_deg_summary <- function(summary_df, out_dir) {
  if (is.null(summary_df) || nrow(summary_df) == 0) return(invisible(NULL))
  
  p1 <- ggplot2::ggplot(summary_df, ggplot2::aes(x = contrast, y = n_sig, fill = analysis)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::geom_text(
      ggplot2::aes(label = n_sig),
      position = ggplot2::position_dodge(width = 0.9),
      vjust = -0.3,
      size = 3
    ) +
    ggplot2::facet_wrap(~ celltype, scales = "free_y") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
    ggplot2::labs(
      title = "DEG counts: pooled vs within-pool vs paired",
      x = "Contrast",
      y = "Number of significant DEGs"
    )
  
  save_plot(p1, file.path(out_dir, "DEG_side_by_side_total"), width = 16, height = 10)
  
  summary_long <- summary_df %>%
    tidyr::pivot_longer(
      cols = c(n_up, n_down),
      names_to = "direction",
      values_to = "n_genes"
    ) %>%
    dplyr::mutate(
      direction = dplyr::case_when(
        direction == "n_up" ~ "Up",
        direction == "n_down" ~ "Down",
        TRUE ~ direction
      )
    )
  
  p2 <- ggplot2::ggplot(summary_long, ggplot2::aes(x = contrast, y = n_genes, fill = direction)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(
      ggplot2::aes(label = n_genes),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3
    ) +
    ggplot2::facet_grid(celltype ~ analysis, scales = "free_y") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
    ggplot2::labs(
      title = "Up- and down-regulated DEG counts side-by-side",
      x = "Contrast",
      y = "Number of genes"
    )
  
  save_plot(p2, file.path(out_dir, "DEG_side_by_side_up_down"), width = 18, height = 12)
}

plot_gsea_summary <- function(summary_df, out_dir) {
  if (is.null(summary_df) || nrow(summary_df) == 0) return(invisible(NULL))
  
  p1 <- ggplot2::ggplot(summary_df, ggplot2::aes(x = contrast, y = n_sig, fill = database)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::geom_text(
      ggplot2::aes(label = n_sig),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3
    ) +
    ggplot2::facet_grid(celltype ~ analysis, scales = "free_y") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
    ggplot2::labs(
      title = "GSEA counts: pooled vs within-pool vs paired",
      x = "Contrast",
      y = "Number of significant pathways"
    )
  
  save_plot(p1, file.path(out_dir, "GSEA_side_by_side_total"), width = 18, height = 12)
  
  p2 <- ggplot2::ggplot(summary_df, ggplot2::aes(x = database, y = n_sig, fill = database)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = n_sig), vjust = -0.3, size = 3) +
    ggplot2::facet_grid(celltype ~ contrast, scales = "free_y") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
    ggplot2::labs(
      title = "Significant pathways by database",
      x = "Database",
      y = "Number of significant pathways"
    )
  
  save_plot(p2, file.path(out_dir, "GSEA_side_by_side_by_database"), width = 20, height = 12)
}

# =============================================================================
# BUILD PSEUDOBULK
# =============================================================================
build_pseudobulk <- function(obj,
                              celltype_col = "Monaco_main_pruned",
                            # celltype_col = "Monaco_main_pruned",
                            #celltype_col = "cluster",
                             assay = "RNA",
                             min_cells_per_pb = 20) {
  log_msg("PB", "Building pseudobulk matrix")
  
  md <- obj@meta.data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")
  
  required_cols <- c("Treatment", "subject", "HTO_pool", celltype_col)
  missing_cols <- setdiff(required_cols, colnames(md))
  if (length(missing_cols) > 0) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
  }
  
  md <- md %>%
    dplyr::mutate(
      CellType = as.character(.data[[celltype_col]]),
      subject_id = as.character(subject),
      Pool = as.character(HTO_pool),
      group = dplyr::case_when(
        Treatment == "HC"        ~ "HC",
        Treatment == "Severe"    ~ "Severe",
        Treatment == "Untreated" ~ "Mild_Untreated",
        Treatment == "Treated"   ~ "Mild_Treated",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(
      !is.na(CellType),
      !is.na(subject_id),
      !is.na(Pool),
      !is.na(group)
    ) %>%
    dplyr::mutate(
      sample_id = paste(subject_id, Pool, group, CellType, sep = "__")
    )
  
  counts_mat <- get_counts_matrix_v5(obj, assay = assay)
  
  common_cells <- intersect(md$cell, colnames(counts_mat))
  md <- md %>% dplyr::filter(cell %in% common_cells)
  counts_mat <- counts_mat[, common_cells, drop = FALSE]
  
  split_cells <- split(md$cell, md$sample_id)
  
  pb_counts_list <- list()
  pb_ncells <- integer(0)
  
  for (sid in names(split_cells)) {
    cells_use <- intersect(split_cells[[sid]], colnames(counts_mat))
    if (length(cells_use) < min_cells_per_pb) next
    pb_counts_list[[sid]] <- Matrix::rowSums(counts_mat[, cells_use, drop = FALSE])
    pb_ncells[sid] <- length(cells_use)
  }
  
  if (length(pb_counts_list) == 0) {
    stop("No pseudobulk samples passed min_cells_per_pb")
  }
  
  pb_counts <- do.call(cbind, pb_counts_list)
  pb_counts <- as.matrix(pb_counts)
  
  sample_meta <- md %>%
    dplyr::select(sample_id, subject_id, Pool, group, CellType) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE)
  
  rownames(sample_meta) <- sample_meta$sample_id
  sample_meta <- sample_meta[colnames(pb_counts), , drop = FALSE]
  sample_meta$n_cells_pb <- pb_ncells[colnames(pb_counts)]
  
  list(counts = pb_counts, meta = sample_meta)
}

# =============================================================================
# LIMMA
# =============================================================================
run_limma_contrast <- function(counts_mat, design, contrast_formula = NULL,
                               coef_name = NULL, out_csv) {
  dge <- edgeR::DGEList(counts = counts_mat)
  keep <- edgeR::filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  if (nrow(dge) == 0) {
    write_skip(sub("\\.csv$", "_SKIPPED.csv", out_csv), "No genes passed filterByExpr")
    return(NULL)
  }
  
  dge <- edgeR::calcNormFactors(dge)
  v <- limma::voom(dge, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  
  if (!is.null(contrast_formula)) {
    cm <- limma::makeContrasts(contrasts = contrast_formula, levels = design)
    fit <- limma::contrasts.fit(fit, cm)
    fit <- limma::eBayes(fit)
    tt <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "P")
  } else {
    fit <- limma::eBayes(fit)
    tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  }
  
  tt$gene <- rownames(tt)
  data.table::fwrite(tt, out_csv)
  tt
}

# =============================================================================
# GSEA
# =============================================================================
run_gsea_for_deg <- function(tt, out_dir, prefix, top_n_plot = 15) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  debug_log <- file.path(out_dir, paste0(prefix, "_GSEA_DEBUG.txt"))
  log_line <- function(...) {
    msg <- paste0(...)
    cat(msg, "\n")
    write(msg, file = debug_log, append = TRUE)
  }
  
  log_line("====================================================")
  log_line("Starting GSEA for: ", prefix)
  
  if (is.null(tt) || nrow(tt) == 0) {
    log_line("SKIP: tt is NULL or empty")
    write_skip(file.path(out_dir, paste0(prefix, "_GSEA_SKIPPED.csv")),
               "DEG table is NULL or empty")
    return(NULL)
  }
  
  if (!all(c("gene", "t") %in% colnames(tt))) {
    log_line("SKIP: missing gene and/or t columns")
    write_skip(file.path(out_dir, paste0(prefix, "_GSEA_SKIPPED.csv")),
               "Missing gene and/or t columns")
    return(NULL)
  }
  
  ranked_sym <- tt %>%
    dplyr::filter(!is.na(gene), gene != "", !is.na(t)) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(t = t[which.max(abs(t))], .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(t))
  
  log_line("Ranked SYMBOL genes: ", nrow(ranked_sym))
  
  if (nrow(ranked_sym) < 20) {
    log_line("SKIP: fewer than 20 ranked genes")
    write_skip(file.path(out_dir, paste0(prefix, "_GSEA_SKIPPED.csv")),
               "Too few ranked genes")
    return(NULL)
  }
  
  gene_list_symbol <- ranked_sym$t
  names(gene_list_symbol) <- ranked_sym$gene
  gene_list_symbol <- sort(gene_list_symbol, decreasing = TRUE)
  
  gene_list_entrez <- tryCatch(
    symbol_to_entrez_ranked(gene_list_symbol),
    error = function(e) {
      log_line("ENTREZ mapping failed: ", e$message)
      NULL
    }
  )
  
  log_line("Mapped ENTREZ genes: ", ifelse(is.null(gene_list_entrez), 0, length(gene_list_entrez)))
  if (!is.null(gene_list_entrez)) {
    log_line("ENTREZ class: ", paste(class(gene_list_entrez), collapse = ","))
    log_line("ENTREZ is.numeric: ", is.numeric(gene_list_entrez))
  }
  
  results <- list()
  
  # 1. GO BP FIRST
  if (!is.null(gene_list_entrez) && length(gene_list_entrez) >= 20) {
    log_line("Running GO BP for: ", prefix)
    
    go_res <- tryCatch({
      clusterProfiler::gseGO(
        geneList = gene_list_entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        minGSSize = GSEA_MIN_GS,
        maxGSSize = GSEA_MAX_GS,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        eps = 0,
        verbose = FALSE
      )
    }, error = function(e) {
      log_line("GO BP FAILED for ", prefix, ": ", e$message)
      NULL
    })
    
    if (!is.null(go_res)) {
      res <- as.data.frame(go_res)
      log_line("GO BP rows: ", nrow(res))
      data.table::fwrite(res, file.path(out_dir, paste0(prefix, "_GO_BP_results.csv")))
      data.table::fwrite(
        res %>% dplyr::filter(!is.na(p.adjust) & p.adjust < GSEA_PADJ_CUTOFF),
        file.path(out_dir, paste0(prefix, "_GO_BP_results_sig.csv"))
      )
      
      tryCatch(
        plot_gsea_bar(
          res %>% dplyr::filter(!is.na(p.adjust) & p.adjust < GSEA_PADJ_CUTOFF),
          paste0(prefix, " - GO BP"),
          file.path(out_dir, paste0(prefix, "_GO_BP_barplot")),
          top_n = top_n_plot
        ),
        error = function(e) {
          log_line("GO BP plot FAILED for ", prefix, ": ", e$message)
        }
      )
      
      results$GO_BP <- res
    } else {
      write_skip(file.path(out_dir, paste0(prefix, "_GO_BP_SKIPPED.csv")),
                 "GO BP returned NULL")
    }
  } else {
    log_line("SKIP GO BP: insufficient ENTREZ mapping")
    write_skip(file.path(out_dir, paste0(prefix, "_GO_BP_SKIPPED.csv")),
               "Insufficient ENTREZ mapping")
  }
  
  # 2. REACTOME SECOND
  if (!is.null(gene_list_entrez) && length(gene_list_entrez) >= 20) {
    log_line("Running Reactome for: ", prefix)
    
    reactome_res <- tryCatch({
      ReactomePA::gsePathway(
        geneList = gene_list_entrez,
        organism = "human",
        minGSSize = GSEA_MIN_GS,
        maxGSSize = GSEA_MAX_GS,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        eps = 0,
        verbose = FALSE
      )
    }, error = function(e) {
      log_line("Reactome FAILED for ", prefix, ": ", e$message)
      NULL
    })
    
    if (!is.null(reactome_res)) {
      res <- as.data.frame(reactome_res)
      log_line("Reactome rows: ", nrow(res))
      data.table::fwrite(res, file.path(out_dir, paste0(prefix, "_Reactome_results.csv")))
      data.table::fwrite(
        res %>% dplyr::filter(!is.na(p.adjust) & p.adjust < GSEA_PADJ_CUTOFF),
        file.path(out_dir, paste0(prefix, "_Reactome_results_sig.csv"))
      )
      
      tryCatch(
        plot_gsea_bar(
          res %>% dplyr::filter(!is.na(p.adjust) & p.adjust < GSEA_PADJ_CUTOFF),
          paste0(prefix, " - Reactome"),
          file.path(out_dir, paste0(prefix, "_Reactome_barplot")),
          top_n = top_n_plot
        ),
        error = function(e) {
          log_line("Reactome plot FAILED for ", prefix, ": ", e$message)
        }
      )
      
      results$Reactome <- res
    } else {
      write_skip(file.path(out_dir, paste0(prefix, "_Reactome_SKIPPED.csv")),
                 "Reactome returned NULL")
    }
  } else {
    log_line("SKIP Reactome: insufficient ENTREZ mapping")
    write_skip(file.path(out_dir, paste0(prefix, "_Reactome_SKIPPED.csv")),
               "Insufficient ENTREZ mapping")
  }
  
  # 3. HALLMARK LAST
  log_line("Running Hallmark for: ", prefix)
  
  hallmark_df <- tryCatch({
    msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
      dplyr::select(gs_name, gene_symbol) %>%
      dplyr::distinct()
  }, error = function(e) {
    log_line("Hallmark gene set load FAILED: ", e$message)
    NULL
  })
  
  if (!is.null(hallmark_df)) {
    hallmark_res <- tryCatch({
      clusterProfiler::GSEA(
        geneList = gene_list_symbol,
        TERM2GENE = hallmark_df,
        minGSSize = GSEA_MIN_GS,
        maxGSSize = GSEA_MAX_GS,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        eps = 0,
        verbose = FALSE
      )
    }, error = function(e) {
      log_line("Hallmark FAILED for ", prefix, ": ", e$message)
      NULL
    })
    
    if (!is.null(hallmark_res)) {
      res <- as.data.frame(hallmark_res)
      log_line("Hallmark rows: ", nrow(res))
      data.table::fwrite(res, file.path(out_dir, paste0(prefix, "_Hallmark_results.csv")))
      data.table::fwrite(
        res %>% dplyr::filter(!is.na(p.adjust) & p.adjust < GSEA_PADJ_CUTOFF),
        file.path(out_dir, paste0(prefix, "_Hallmark_results_sig.csv"))
      )
      
      tryCatch(
        plot_gsea_bar(
          res %>% dplyr::filter(!is.na(p.adjust) & p.adjust < GSEA_PADJ_CUTOFF),
          paste0(prefix, " - Hallmark"),
          file.path(out_dir, paste0(prefix, "_Hallmark_barplot")),
          top_n = top_n_plot
        ),
        error = function(e) {
          log_line("Hallmark plot FAILED for ", prefix, ": ", e$message)
        }
      )
      
      results$Hallmark <- res
    } else {
      write_skip(file.path(out_dir, paste0(prefix, "_Hallmark_SKIPPED.csv")),
                 "Hallmark returned NULL")
    }
  } else {
    write_skip(file.path(out_dir, paste0(prefix, "_Hallmark_SKIPPED.csv")),
               "Hallmark gene sets failed to load")
  }
  
  log_line("Finished GSEA for: ", prefix)
  log_line("Databases returned: ", paste(names(results), collapse = ", "))
  
  if (length(results) == 0) {
    log_line("No databases returned results")
    return(NULL)
  }
  
  results
}

summarize_gsea_results <- function(gsea_res, celltype, analysis, contrast) {
  if (is.null(gsea_res)) return(NULL)
  
  out <- list()
  for (db in names(gsea_res)) {
    df <- gsea_res[[db]]
    if (is.null(df) || nrow(df) == 0) {
      out[[length(out) + 1]] <- data.frame(
        celltype = celltype,
        analysis = analysis,
        contrast = contrast,
        database = db,
        n_sig = 0,
        stringsAsFactors = FALSE
      )
    } else {
      out[[length(out) + 1]] <- data.frame(
        celltype = celltype,
        analysis = analysis,
        contrast = contrast,
        database = db,
        n_sig = sum(df$p.adjust < GSEA_PADJ_CUTOFF, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  dplyr::bind_rows(out)
}

# =============================================================================
# ONE ANALYSIS BLOCK
# =============================================================================
run_unpaired_analysis <- function(ct, ct_counts, ct_meta, out_dir, gsea_dir,
                                  analysis_label, groups_keep, contrast_name,
                                  contrast_formula, pool_subset = NULL) {
  meta_use <- ct_meta %>%
    dplyr::filter(group %in% groups_keep)
  
  if (!is.null(pool_subset)) {
    meta_use <- meta_use %>%
      dplyr::filter(Pool %in% pool_subset)
  }
  
  grp_n <- table(meta_use$group)
  needed <- unique(groups_keep)
  
  out_csv <- file.path(out_dir, paste0("DEG_", safe_name(ct), "_", contrast_name, ".csv"))
  
  if (!all(needed %in% names(grp_n)) || any(grp_n[needed] < MIN_SAMPLES_GROUP)) {
    write_skip(out_csv, paste0("Insufficient samples for ", contrast_name))
    return(list(deg = NULL, gsea = NULL, deg_summary = NULL, gsea_summary = NULL))
  }
  
  meta_use$group <- factor(meta_use$group, levels = groups_keep)
  design <- model.matrix(~ 0 + group, data = meta_use)
  colnames(design) <- make.names(colnames(design))
  
  tt <- tryCatch(
    run_limma_contrast(
      counts_mat = ct_counts[, rownames(meta_use), drop = FALSE],
      design = design,
      contrast_formula = contrast_formula,
      out_csv = out_csv
    ),
    error = function(e) {
      write_skip(out_csv, paste0("limma failed: ", e$message))
      NULL
    }
  )
  
  if (is.null(tt)) {
    return(list(deg = NULL, gsea = NULL, deg_summary = NULL, gsea_summary = NULL))
  }
  
  deg_summary <- data.frame(
    celltype = ct,
    analysis = analysis_label,
    contrast = contrast_name,
    n_up = sum(tt$adj.P.Val < DEG_FDR_CUTOFF & tt$logFC >  DEG_LOGFC_CUTOFF, na.rm = TRUE),
    n_down = sum(tt$adj.P.Val < DEG_FDR_CUTOFF & tt$logFC < -DEG_LOGFC_CUTOFF, na.rm = TRUE),
    n_sig = sum(tt$adj.P.Val < DEG_FDR_CUTOFF & abs(tt$logFC) >= DEG_LOGFC_CUTOFF, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  gsea_res <- run_gsea_for_deg(tt, gsea_dir, paste0(safe_name(ct), "_", contrast_name), top_n_plot = GSEA_TOP_N_PLOT)
  gsea_summary <- summarize_gsea_results(gsea_res, ct, analysis_label, contrast_name)
  
  list(deg = tt, gsea = gsea_res, deg_summary = deg_summary, gsea_summary = gsea_summary)
}

# FIXED: handles missing Mild_Untreated / Mild_Treated columns safely
run_paired_analysis <- function(ct, ct_counts, ct_meta, out_dir, gsea_dir) {
  out_csv <- file.path(out_dir, paste0("DEG_", safe_name(ct), "_MildTreated_vs_MildUntreated_paired.csv"))
  
  meta_use <- ct_meta %>%
    dplyr::filter(group %in% c("Mild_Untreated", "Mild_Treated"))
  
  if (nrow(meta_use) == 0) {
    write_skip(out_csv, "No mild paired samples for this cell type")
    return(list(deg = NULL, gsea = NULL, deg_summary = NULL, gsea_summary = NULL))
  }
  
  pair_tab <- meta_use %>%
    dplyr::count(subject_id, group) %>%
    tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0)
  
  if (!("Mild_Untreated" %in% colnames(pair_tab))) pair_tab$Mild_Untreated <- 0
  if (!("Mild_Treated" %in% colnames(pair_tab))) pair_tab$Mild_Treated <- 0
  
  pair_subjects <- pair_tab %>%
    dplyr::filter(Mild_Untreated > 0, Mild_Treated > 0) %>%
    dplyr::pull(subject_id)
  
  meta_use <- meta_use %>%
    dplyr::filter(subject_id %in% pair_subjects) %>%
    dplyr::arrange(subject_id, group)
  
  if (length(unique(meta_use$subject_id)) < 2) {
    write_skip(out_csv, "Fewer than 2 complete paired subjects")
    return(list(deg = NULL, gsea = NULL, deg_summary = NULL, gsea_summary = NULL))
  }
  
  meta_use$group <- factor(meta_use$group, levels = c("Mild_Untreated", "Mild_Treated"))
  meta_use$subject_id <- factor(meta_use$subject_id)
  
  design <- model.matrix(~ subject_id + group, data = meta_use)
  coef_name <- grep("^group", colnames(design), value = TRUE)[1]
  
  tt <- tryCatch(
    run_limma_contrast(
      counts_mat = ct_counts[, rownames(meta_use), drop = FALSE],
      design = design,
      coef_name = coef_name,
      out_csv = out_csv
    ),
    error = function(e) {
      write_skip(out_csv, paste0("limma failed: ", e$message))
      NULL
    }
  )
  
  if (is.null(tt)) {
    return(list(deg = NULL, gsea = NULL, deg_summary = NULL, gsea_summary = NULL))
  }
  
  contrast_name <- "MildTreated_vs_MildUntreated_paired"
  
  deg_summary <- data.frame(
    celltype = ct,
    analysis = "paired_cross_pool",
    contrast = contrast_name,
    n_up = sum(tt$adj.P.Val < DEG_FDR_CUTOFF & tt$logFC >  DEG_LOGFC_CUTOFF, na.rm = TRUE),
    n_down = sum(tt$adj.P.Val < DEG_FDR_CUTOFF & tt$logFC < -DEG_LOGFC_CUTOFF, na.rm = TRUE),
    n_sig = sum(tt$adj.P.Val < DEG_FDR_CUTOFF & abs(tt$logFC) >= DEG_LOGFC_CUTOFF, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  gsea_res <- run_gsea_for_deg(tt, gsea_dir, paste0(safe_name(ct), "_", contrast_name), top_n_plot = GSEA_TOP_N_PLOT)
  gsea_summary <- summarize_gsea_results(gsea_res, ct, "paired_cross_pool", contrast_name)
  
  list(deg = tt, gsea = gsea_res, deg_summary = deg_summary, gsea_summary = gsea_summary)
}

# =============================================================================
# MAIN DRIVER
# Added:
#   start_from_celltype = NULL   -> run all cell types
#   start_from_celltype = "Monocytes" -> start at Monocytes and continue onward
# =============================================================================
run_deg_gsea_pooled_and_within_pool <- function(obj,
                                                out_dir = "CLUSTER_deg_gsea_pooled_vs_within_pool",
                                              #  celltype_col = "Monaco_main_pruned",
                                              # celltype_col = "cluster",
                                              celltype_col = "Monaco_main_pruned",
                                                assay = "RNA",
                                                start_from_celltype = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  pooled_deg_dir   <- file.path(out_dir, "pooled", "deg_by_celltype")
  pooled_gsea_dir  <- file.path(out_dir, "pooled", "gsea_by_celltype")
  within_deg_dir   <- file.path(out_dir, "within_pool", "deg_by_celltype")
  within_gsea_dir  <- file.path(out_dir, "within_pool", "gsea_by_celltype")
  paired_deg_dir   <- file.path(out_dir, "paired_cross_pool", "deg_by_celltype")
  paired_gsea_dir  <- file.path(out_dir, "paired_cross_pool", "gsea_by_celltype")
  
  dir.create(pooled_deg_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pooled_gsea_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(within_deg_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(within_gsea_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paired_deg_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(paired_gsea_dir, recursive = TRUE, showWarnings = FALSE)
  
  log_msg("DEG", "Section 4: Pseudobulk DEG + GSEA side-by-side")
  
  pb <- build_pseudobulk(
    obj = obj,
    celltype_col = celltype_col,
    assay = assay,
    min_cells_per_pb = MIN_CELLS_PER_PB
  )
  
  pb_counts <- pb$counts
  sample_meta <- pb$meta
  
  data.table::fwrite(sample_meta, file.path(out_dir, "pseudobulk_sample_metadata.csv"))
  
  all_deg_summaries <- list()
  all_gsea_summaries <- list()
  
  celltypes <- sort(unique(sample_meta$CellType))
  
  if (!is.null(start_from_celltype)) {
    if (!(start_from_celltype %in% celltypes)) {
      stop("start_from_celltype not found in cell types: ", start_from_celltype)
    }
    start_idx <- match(start_from_celltype, celltypes)
    celltypes <- celltypes[start_idx:length(celltypes)]
    log_msg("DEG", paste("Resuming from cell type:", start_from_celltype))
  }
  
  for (ct in celltypes) {
    log_msg("DEG", paste("  Processing", ct))
    ct_safe <- safe_name(ct)
    
    ct_meta <- sample_meta %>%
      dplyr::filter(CellType == ct, n_cells_pb >= MIN_CELLS_PER_PB)
    
    if (nrow(ct_meta) == 0) next
    
    ct_counts <- pb_counts[, rownames(ct_meta), drop = FALSE]
    
    ct_pooled_deg_dir  <- file.path(pooled_deg_dir, ct_safe)
    ct_pooled_gsea_dir <- file.path(pooled_gsea_dir, ct_safe)
    ct_within_deg_dir  <- file.path(within_deg_dir, ct_safe)
    ct_within_gsea_dir <- file.path(within_gsea_dir, ct_safe)
    ct_paired_deg_dir  <- file.path(paired_deg_dir, ct_safe)
    ct_paired_gsea_dir <- file.path(paired_gsea_dir, ct_safe)
    
    dir.create(ct_pooled_deg_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(ct_pooled_gsea_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(ct_within_deg_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(ct_within_gsea_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(ct_paired_deg_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(ct_paired_gsea_dir, recursive = TRUE, showWarnings = FALSE)
    
    # POOLED
    res1 <- run_unpaired_analysis(
      ct = ct,
      ct_counts = ct_counts,
      ct_meta = ct_meta,
      out_dir = ct_pooled_deg_dir,
      gsea_dir = ct_pooled_gsea_dir,
      analysis_label = "pooled",
      groups_keep = c("Mild_Untreated", "HC"),
      contrast_name = "HC_vs_MildUntreated_pooled",
      contrast_formula = "groupHC - groupMild_Untreated",
      pool_subset = NULL
    )
    if (!is.null(res1$deg_summary))  all_deg_summaries[[length(all_deg_summaries) + 1]] <- res1$deg_summary
    if (!is.null(res1$gsea_summary)) all_gsea_summaries[[length(all_gsea_summaries) + 1]] <- res1$gsea_summary
    
    res2 <- run_unpaired_analysis(
      ct = ct,
      ct_counts = ct_counts,
      ct_meta = ct_meta,
      out_dir = ct_pooled_deg_dir,
      gsea_dir = ct_pooled_gsea_dir,
      analysis_label = "pooled",
      groups_keep = c("Mild_Untreated", "Severe"),
      contrast_name = "Severe_vs_MildUntreated_pooled",
      contrast_formula = "groupSevere - groupMild_Untreated",
      pool_subset = NULL
    )
    if (!is.null(res2$deg_summary))  all_deg_summaries[[length(all_deg_summaries) + 1]] <- res2$deg_summary
    if (!is.null(res2$gsea_summary)) all_gsea_summaries[[length(all_gsea_summaries) + 1]] <- res2$gsea_summary
    
    res3 <- run_unpaired_analysis(
      ct = ct,
      ct_counts = ct_counts,
      ct_meta = ct_meta,
      out_dir = ct_pooled_deg_dir,
      gsea_dir = ct_pooled_gsea_dir,
      analysis_label = "pooled",
      groups_keep = c("HC", "Severe"),
      contrast_name = "Severe_vs_HC_pooled",
      contrast_formula = "groupSevere - groupHC",
      pool_subset = NULL
    )
    if (!is.null(res3$deg_summary))  all_deg_summaries[[length(all_deg_summaries) + 1]] <- res3$deg_summary
    if (!is.null(res3$gsea_summary)) all_gsea_summaries[[length(all_gsea_summaries) + 1]] <- res3$gsea_summary
    
    # WITHIN POOL1
    res4 <- run_unpaired_analysis(
      ct = ct,
      ct_counts = ct_counts,
      ct_meta = ct_meta,
      out_dir = ct_within_deg_dir,
      gsea_dir = ct_within_gsea_dir,
      analysis_label = "within_pool1",
      groups_keep = c("Mild_Untreated", "HC"),
      contrast_name = "HC_vs_MildUntreated_Pool1",
      contrast_formula = "groupHC - groupMild_Untreated",
      pool_subset = "Pool1"
    )
    if (!is.null(res4$deg_summary))  all_deg_summaries[[length(all_deg_summaries) + 1]] <- res4$deg_summary
    if (!is.null(res4$gsea_summary)) all_gsea_summaries[[length(all_gsea_summaries) + 1]] <- res4$gsea_summary
    
    res5 <- run_unpaired_analysis(
      ct = ct,
      ct_counts = ct_counts,
      ct_meta = ct_meta,
      out_dir = ct_within_deg_dir,
      gsea_dir = ct_within_gsea_dir,
      analysis_label = "within_pool1",
      groups_keep = c("Mild_Untreated", "Severe"),
      contrast_name = "Severe_vs_MildUntreated_Pool1",
      contrast_formula = "groupSevere - groupMild_Untreated",
      pool_subset = "Pool1"
    )
    if (!is.null(res5$deg_summary))  all_deg_summaries[[length(all_deg_summaries) + 1]] <- res5$deg_summary
    if (!is.null(res5$gsea_summary)) all_gsea_summaries[[length(all_gsea_summaries) + 1]] <- res5$gsea_summary
    
    res6 <- run_unpaired_analysis(
      ct = ct,
      ct_counts = ct_counts,
      ct_meta = ct_meta,
      out_dir = ct_within_deg_dir,
      gsea_dir = ct_within_gsea_dir,
      analysis_label = "within_pool1",
      groups_keep = c("HC", "Severe"),
      contrast_name = "Severe_vs_HC_Pool1",
      contrast_formula = "groupSevere - groupHC",
      pool_subset = "Pool1"
    )
    if (!is.null(res6$deg_summary))  all_deg_summaries[[length(all_deg_summaries) + 1]] <- res6$deg_summary
    if (!is.null(res6$gsea_summary)) all_gsea_summaries[[length(all_gsea_summaries) + 1]] <- res6$gsea_summary
    
    # WITHIN POOL2
    res7 <- run_unpaired_analysis(
      ct = ct,
      ct_counts = ct_counts,
      ct_meta = ct_meta,
      out_dir = ct_within_deg_dir,
      gsea_dir = ct_within_gsea_dir,
      analysis_label = "within_pool2",
      groups_keep = c("HC", "Severe"),
      contrast_name = "Severe_vs_HC_Pool2",
      contrast_formula = "groupSevere - groupHC",
      pool_subset = "Pool2"
    )
    if (!is.null(res7$deg_summary))  all_deg_summaries[[length(all_deg_summaries) + 1]] <- res7$deg_summary
    if (!is.null(res7$gsea_summary)) all_gsea_summaries[[length(all_gsea_summaries) + 1]] <- res7$gsea_summary
    
    # PAIRED CROSS-POOL
    res8 <- tryCatch(
      run_paired_analysis(
        ct = ct,
        ct_counts = ct_counts,
        ct_meta = ct_meta,
        out_dir = ct_paired_deg_dir,
        gsea_dir = ct_paired_gsea_dir
      ),
      error = function(e) {
        message("Paired analysis failed for ", ct, ": ", e$message)
        NULL
      }
    )
    if (!is.null(res8)) {
      if (!is.null(res8$deg_summary))  all_deg_summaries[[length(all_deg_summaries) + 1]] <- res8$deg_summary
      if (!is.null(res8$gsea_summary)) all_gsea_summaries[[length(all_gsea_summaries) + 1]] <- res8$gsea_summary
    }
  }
  
  deg_summary <- dplyr::bind_rows(all_deg_summaries)
  gsea_summary <- dplyr::bind_rows(all_gsea_summaries)
  
  if (nrow(deg_summary) > 0) {
    data.table::fwrite(deg_summary, file.path(out_dir, "DEG_side_by_side_summary.csv"))
    plot_deg_summary(deg_summary, out_dir)
  }
  
  if (nrow(gsea_summary) > 0) {
    data.table::fwrite(gsea_summary, file.path(out_dir, "GSEA_side_by_side_summary.csv"))
    plot_gsea_summary(gsea_summary, out_dir)
  }
  
  invisible(list(
    pseudobulk = pb,
    deg_summary = deg_summary,
    gsea_summary = gsea_summary
  ))
}

# Example full run:
# res_side_by_side <- run_deg_gsea_pooled_and_within_pool(
#   obj = obj_filtered,
#   out_dir = "deg_gsea_pooled_vs_within_pool",
#   celltype_col = "Monaco_main_pruned",
#   assay = "RNA"
# )

# Example resume from Monocytes onward:
# res_side_by_side <- run_deg_gsea_pooled_and_within_pool(
#   obj = obj_filtered,
#   out_dir = "deg_gsea_pooled_vs_within_pool_resume",
#   celltype_col = "Monaco_main_pruned",
#   assay = "RNA",
#   start_from_celltype = "Monocytes"
# )


#res_side_by_side <- run_deg_gsea_pooled_and_within_pool(
#  obj = obj_filtered,
#  out_dir = "deg_gsea_pooled_vs_within_pool_resume",
#  celltype_col = "Monaco_main_pruned",
#  assay = "RNA",
#  start_from_celltype = "Neutrophils"
#)

#obj_filtered<-readRDS("./PBMC_RDS/object_final_annotated.rds")

# Example full run:
# res_side_by_side <- run_deg_gsea_pooled_and_within_pool(
#   obj = obj_filtered,
#   out_dir = "deg_gsea_pooled_vs_within_pool",
#   celltype_col = "cluster",
#   assay = "RNA"
# )

##"Monaco_fine_pruned"

 ##Example full run:
 res_side_by_side <- run_deg_gsea_pooled_and_within_pool(
   obj = obj_filtered,
   out_dir = "deg_gsea_pooled_vs_within_pool_Monaco_main_pruned",
   celltype_col = "Monaco_main_pruned",
   assay = "RNA"
 )
