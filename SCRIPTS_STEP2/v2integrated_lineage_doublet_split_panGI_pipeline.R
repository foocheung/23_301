#!/usr/bin/env Rscript
## ============================================================================
## step5_6_merged.R  (UPDATED: supports PanGI_L1 split)
##
## Purpose:
##   Minimal-change merge of step5.R and step6.R into a single script.
##   1) Loads a Seurat object, maps clusters -> lineage, preserves global UMAP,
##      generates global plots.
##   2) Splits by lineage, reclusters each subset, saves per-lineage RDS + plots.
##   3) (Optional, default TRUE) Runs doublet removal & re-clustering on each
##      per-lineage object using the logic from step6.R.
##
## New:
##   --split_by csv|pangi
##     csv   = original behavior using annot_csv
##     pangi = set lineage_call from obj$PanGI_L1:
##             Haemopoietic = {T and NK cells, B and B plasma}
##             Non-Haemopoietic = {Myeloid, Epithelial, Endothelial, Mesenchymal, Neural}
##             Unknown if none of the above
## ============================================================================
## Useus this paper https://www.nature.com/articles/s41586-024-07571-1 gutATLAS
suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(dplyr)
  library(Matrix)
  library(data.table)
  library(stringr)
})

## -------------------------- CLI ARGUMENTS -----------------------------------
opt_list <- list(
  make_option("--in_rds", type = "character",
              default = "./MERGED_PIPELINE_OUTPUT_V5/object_final_annotated_panGI.rds",
              help = "Input Seurat RDS [default: %default]"),
  make_option("--annot_csv", type = "character",
              default = "./MERGED_PIPELINE_OUTPUT_V5/cluster_annotations_comprehensive.csv",
              help = "Annotation CSV with columns: cluster,lineage_call [default: %default]"),
  make_option("--outdir_step5", type = "character", default = "LINEAGE_DOUBLET_ANALYSIS_PANGI_V9",
              help = "Output directory for Step 5 products [default: %default]"),
  make_option("--assay", type = "character", default = "RNA",
              help = "Default assay to use [default: %default]"),
  make_option("--use_sct", action = "store_true", default = FALSE,
              help = "Use SCTransform for subset reclustering [default: %default]"),
  make_option("--n_pcs", type = "integer", default = 50,
              help = "Number of PCs [default: %default]"),
  make_option("--umap_dims", type = "character", default = "1:30",
              help = "UMAP dims (R-style, e.g. '1:30') [default: %default]"),
  make_option("--resolution5", type = "double", default = 0.6,
              help = "Resolution for subset FindClusters (Step 5) [default: %default]"),
  ## Step6-style options
  make_option("--run_step6", action = "store_true", default = TRUE,
              help = "Run doublet removal (Step 6) on all per-lineage RDS [default: %default]"),
  make_option("--outdir_step6", type = "character", default = "LINEAGE_DOUBLET_ANALYSIS_PANGI_V9",
              help = "Output directory for Step 6 products [default: %default]"),
  make_option("--dims6", type = "integer", default = 30,
              help = "Number of PCs for Step 6 clustering [default: %default]"),
  make_option("--resolution6", type = "character", default = "0.3,0.5,0.8,1.0",
              help = "Comma-separated clustering resolutions for Step 6 [default: %default]"),
  make_option("--cluster_doublet_threshold", type = "double", default = 70,
              help = "Remove entire clusters with >= this %% doublets [default: %default]"),
  ## NEW: split mode
  make_option("--split_by", type = "character", default = "csv",
              help = "Split mode: 'csv' (use annot_csv) or 'pangi' (use PanGI_L1 to Haemopoietic/Non-Haemopoietic) [default: %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

## --------------------------- SHARED HELPERS ---------------------------------
section <- function(title) cat(sprintf("\n========== %s ==========\n", title))
logmsg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

detect_cluster_col <- function(obj) {
  if ("seurat_clusters" %in% colnames(obj@meta.data)) return("seurat_clusters")
  cand <- grep("_snn_res", colnames(obj@meta.data), value = TRUE)
  if (length(cand)) return(cand[1])
  stop("No cluster column found in object.")
}

# Seurat v5-friendly counts retrieval with slot fallback
get_counts_matrix <- function(obj, assay) {
  DefaultAssay(obj) <- if (assay %in% names(Assays(obj))) assay else DefaultAssay(obj)
  tryCatch(
    GetAssayData(obj, assay = DefaultAssay(obj), layer = "counts"),
    error = function(e) GetAssayData(obj, assay = DefaultAssay(obj), slot = "counts")
  )
}

ensure_percent_mt <- function(obj) {
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    mt_genes <- grep("^MT-|^Mt-|^mt-", rownames(obj), value = TRUE)
    if (length(mt_genes) > 0) {
      obj <- PercentageFeatureSet(obj, pattern = "^MT-|^Mt-|^mt-", col.name = "percent.mt")
      logmsg("Computed percent.mt (%d MT genes)", length(mt_genes))
    } else {
      obj$percent.mt <- 0
      logmsg("No MT genes; percent.mt set to 0")
    }
  }
  obj
}

ensure_epi_score <- function(obj) {
  if (!"EpithelialScore" %in% colnames(obj@meta.data)) {
    epi_genes <- intersect(c("EPCAM","KRT8","KRT18","KRT19","KRT7"), rownames(obj))
    if (length(epi_genes) > 0) {
      obj <- AddModuleScore(obj, features=list(epi_genes), name="EpithelialScore")
      obj$EpithelialScore <- obj$EpithelialScore1; obj$EpithelialScore1 <- NULL
      logmsg("Computed EpithelialScore (%d genes)", length(epi_genes))
    } else {
      obj$EpithelialScore <- 0
      logmsg("No epithelial marker genes found; EpithelialScore=0")
    }
  }
  obj
}

ensure_cd45_expr <- function(obj) {
  if (!"CD45_expr" %in% colnames(obj@meta.data)) {
    cd45_gene <- intersect(c("PTPRC","Ptprc","ptprc"), rownames(obj))[1]
    if (is.na(cd45_gene)) {
      obj$CD45_expr <- 0
      attr(obj, "cd45_gene") <- NA_character_
      logmsg("CD45/PTPRC not found; CD45_expr set to 0")
    } else {
      obj$CD45_expr <- FetchData(obj, cd45_gene)[,1]
      attr(obj, "cd45_gene") <- cd45_gene
      logmsg("CD45_expr from %s", cd45_gene)
    }
  } else {
    if (is.null(attr(obj, "cd45_gene"))) {
      cand <- intersect(c("PTPRC","Ptprc","ptprc"), rownames(obj))
      attr(obj, "cd45_gene") <- if (length(cand)) cand[1] else NA_character_
    }
  }
  obj
}

ensure_ncount_nfeature <- function(obj, assay = "RNA") {
  counts <- get_counts_matrix(obj, assay)
  if (!"nFeature_RNA" %in% colnames(obj@meta.data)) {
    obj$nFeature_RNA <- Matrix::colSums(counts > 0)
    logmsg("Computed nFeature_RNA")
  }
  if (!"nCount_RNA" %in% colnames(obj@meta.data)) {
    obj$nCount_RNA <- Matrix::colSums(counts)
    logmsg("Computed nCount_RNA")
  }
  obj
}

ensure_umap <- function(obj, assay = "RNA", n_pcs = 50, umap_dims = 1:30) {
  if (!"umap" %in% names(obj@reductions)) {
    logmsg("No UMAP found; computing quick UMAP on %s", DefaultAssay(obj))
    DefaultAssay(obj) <- if (assay %in% names(Assays(obj))) assay else DefaultAssay(obj)
    rn <- obj[[DefaultAssay(obj)]]
    has_data <- tryCatch("data" %in% SeuratObject::Layers(rn), error=function(e) FALSE)
    if (!has_data) obj <- NormalizeData(obj, verbose=FALSE)
    obj <- FindVariableFeatures(obj, nfeatures=3000, verbose=FALSE)
    obj <- ScaleData(obj, verbose=FALSE)
    obj <- RunPCA(obj, npcs=n_pcs, verbose=FALSE)
    obj <- RunUMAP(obj, dims=umap_dims, verbose=FALSE)
  }
  obj
}

recluster_subset <- function(o, assay = "RNA", use_sct = FALSE,
                             n_pcs = 50, umap_dims = 1:30, resolution = 0.6) {
  DefaultAssay(o) <- if (assay %in% names(Assays(o))) assay else DefaultAssay(o)
  if (use_sct) {
    o <- SCTransform(o, vst.flavor="v2", vars.to.regress=c("percent.mt"), verbose=FALSE)
    o <- RunPCA(o, npcs=n_pcs, verbose=FALSE)
  } else {
    rn <- o[[DefaultAssay(o)]]
    has_data <- tryCatch("data" %in% SeuratObject::Layers(rn), error=function(e) FALSE)
    if (!has_data) o <- NormalizeData(o, verbose=FALSE)
    o <- FindVariableFeatures(o, nfeatures=3000, verbose=FALSE)
    o <- ScaleData(o, verbose=FALSE)
    o <- RunPCA(o, npcs=n_pcs, verbose=FALSE)
  }
  o <- RunUMAP(o, dims=umap_dims, verbose=FALSE)
  o <- FindNeighbors(o, dims=umap_dims, verbose=FALSE)
  o <- FindClusters(o, resolution=resolution, verbose=FALSE)
  o
}

feature_umaps <- function(o, prefix, outdir,
                          UMAP_POINT_SIZE=0.4, PLOT_WIDTH=10, PLOT_HEIGHT=8, DPI=300) {
  cd45_gene <- attr(o, "cd45_gene")
  if (is.null(cd45_gene)) {
    cand <- intersect(c("PTPRC","Ptprc","ptprc"), rownames(o))
    cd45_gene <- if (length(cand)) cand[1] else NA_character_
  }
  
  feats <- c("CD45_expr","EpithelialScore","percent.mt","nCount_RNA","nFeature_RNA")
  nmaps <- list(
    CD45_expr      = sprintf("UMAP_%s_feature_CD45.png", prefix),
    EpithelialScore= sprintf("UMAP_%s_feature_EpithelialScore.png", prefix),
    `percent.mt`   = sprintf("UMAP_%s_feature_percent_mt.png", prefix),
    nCount_RNA     = sprintf("UMAP_%s_feature_nCount_RNA.png", prefix),
    nFeature_RNA   = sprintf("UMAP_%s_feature_nFeature_RNA.png", prefix)
  )
  
  for (f in feats) {
    if (!f %in% colnames(o@meta.data)) next
    title_txt <- if (identical(f, "CD45_expr") && !is.null(cd45_gene) && !is.na(cd45_gene)) {
      sprintf("CD45 (%s)", cd45_gene)
    } else f
    
    p <- FeaturePlot(
      o, features = f, reduction = "umap",
      pt.size = UMAP_POINT_SIZE, order = TRUE, raster = TRUE
    ) + ggtitle(title_txt) + theme_minimal()
    
    ggsave(file.path(outdir, nmaps[[f]]), p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = DPI)
  }
}

violin_core <- function(o, prefix, outdir,
                        UMAP_LABEL_SIZE=4, PLOT_WIDTH=10, PLOT_HEIGHT=8, DPI=300) {
  feats <- c("CD45_expr","EpithelialScore","percent.mt","nCount_RNA","nFeature_RNA")
  feats <- feats[feats %in% colnames(o@meta.data)]
  if (length(feats) == 0) return(invisible(NULL))
  Idents(o) <- detect_cluster_col(o)
  for (f in feats) {
    p <- VlnPlot(o, features=f, pt.size=0.05) +
      ggtitle(sprintf("%s by cluster", f)) + theme_minimal()
    ggsave(file.path(outdir, sprintf("Vln_%s_%s.png", prefix, f)),
           p, width=PLOT_WIDTH, height=PLOT_HEIGHT, dpi=DPI)
  }
}

## ----------------------------- STEP 5 ----------------------------------------
run_step5 <- function(in_rds, annot_csv, out_dir, assay="RNA",
                      use_sct=FALSE, n_pcs=50, umap_dims=1:30, resolution=0.6) {
  
  UMAP_POINT_SIZE <- 0.4
  UMAP_LABEL_SIZE <- 4
  PLOT_WIDTH  <- 10
  PLOT_HEIGHT <- 8
  DPI <- 300
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "global"), showWarnings = FALSE)
  dir.create(file.path(out_dir, "subsets"), showWarnings = FALSE)
  
  section("LOAD GLOBAL OBJECT")
  stopifnot(file.exists(in_rds))
  obj <- readRDS(in_rds)
  stopifnot(inherits(obj, "Seurat"))
  logmsg("Loaded object: %d cells, %d features", ncol(obj), nrow(obj))
  
  obj <- ensure_percent_mt(obj)
  obj <- ensure_ncount_nfeature(obj, assay = assay)
  obj <- ensure_epi_score(obj)
  obj <- ensure_cd45_expr(obj)
  obj <- ensure_umap(obj, assay = assay, n_pcs = n_pcs, umap_dims = umap_dims)
  
  if (!"umap_global" %in% names(obj@reductions)) {
    obj[["umap_global"]] <- obj[["umap"]]
    try({ Key(obj[["umap_global"]]) <- "umapg_" }, silent = TRUE)
    logmsg("Stored original UMAP as 'umap_global'")
  }
  
  clust_col <- detect_cluster_col(obj)
  Idents(obj) <- clust_col
  cell_clusters <- as.character(obj@meta.data[[clust_col]])
  
  ## ------------------ CHANGED: lineage_call from CSV OR PanGI_L1 -------------
  section("SET lineage_call (CSV vs PanGI_L1)")
  if (identical(opt$split_by, "pangi")) {
    if (!"PanGI_L1" %in% colnames(obj@meta.data)) {
      stop("split_by='pangi' requested but obj$PanGI_L1 not found. Add Pan-GI labels first.")
    }
    ##"Epithelial"     "T and NK cells" "B and B plasma" "Myeloid"        "Endothelial"    "Mesenchymal"    "Neural" 
    Haemopoietic_vals <- c("T and NK cells","B and B plasma", "Myeloid")
    non_Haemopoietic_vals <- c("Epithelial","Endothelial","Mesenchymal","Neural", "Neural")
    obj$lineage_call <- ifelse(
      obj$PanGI_L1 %in% Haemopoietic_vals, "Haemopoietic",
      ifelse(obj$PanGI_L1 %in% non_Haemopoietic_vals, "Non-Haemopoietic", "Unknown")
    )
    logmsg("Using PanGI_L1 to define lineage_call (Haemopoietic vs Non-Haemopoietic)")
  } else {
    stopifnot(file.exists(annot_csv))
    annot_df <- fread(annot_csv)
    annot_df$cluster <- as.character(annot_df$cluster)
    lineage_map <- setNames(annot_df$lineage_call, annot_df$cluster)
    obj$lineage_call <- unname(lineage_map[cell_clusters])
    obj$lineage_call[is.na(obj$lineage_call)] <- "Unknown"
    logmsg("Using annot_csv to define lineage_call")
  }
  ## ---------------------------------------------------------------------------
  
  lineages <- sort(unique(obj$lineage_call))
  logmsg("Lineages found: %s", paste(lineages, collapse=", "))
  
  section("GLOBAL PLOTS")
  lineage_colors <- c("Haemopoietic"="#E41A1C","Non-Haemopoietic"="#377EB8","Mixed"="#4DAF4A","Unknown"="#999999")
  
  p_global_clust <- DimPlot(obj, reduction="umap_global", group.by=clust_col,
                            label=TRUE, label.size=UMAP_LABEL_SIZE,
                            pt.size=UMAP_POINT_SIZE, repel=TRUE, raster=TRUE) +
    ggtitle("Global UMAP (original) — Clusters") + theme_minimal()
  ggsave(file.path(out_dir, "global", "UMAP_global_clusters.png"),
         p_global_clust, width=PLOT_WIDTH, height=PLOT_HEIGHT, dpi=DPI)
  
  p_global_lineage <- DimPlot(obj, reduction="umap_global", group.by="lineage_call",
                              cols=lineage_colors, pt.size=UMAP_POINT_SIZE, raster=TRUE) +
    ggtitle("Global UMAP (original) — Lineage") + theme_minimal()
  ggsave(file.path(out_dir, "global", "UMAP_global_lineage.png"),
         p_global_lineage, width=PLOT_WIDTH, height=PLOT_HEIGHT, dpi=DPI)
  
  feature_umaps(obj, "global", file.path(out_dir, "global"),
                UMAP_POINT_SIZE = UMAP_POINT_SIZE, PLOT_WIDTH=PLOT_WIDTH, PLOT_HEIGHT=PLOT_HEIGHT, DPI=DPI)
  
  section("SPLIT BY LINEAGE & RECLUSTER")
  subset_paths <- list()
  subset_umaps  <- list()
  
  for (lin in lineages) {
    out_subdir <- file.path(out_dir, "subsets", paste0("lineage_", gsub("[^A-Za-z0-9]+","_", lin)))
    dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
    
    logmsg("Subsetting lineage: %s", lin)
    obj_sub <- subset(obj, subset = lineage_call == lin)
    
    obj_sub <- ensure_percent_mt(obj_sub)
    obj_sub <- ensure_ncount_nfeature(obj_sub, assay = assay)
    obj_sub <- ensure_epi_score(obj_sub)
    obj_sub <- ensure_cd45_expr(obj_sub)
    
    obj_sub <- recluster_subset(obj_sub, assay=assay, use_sct=use_sct,
                                n_pcs=n_pcs, umap_dims=umap_dims, resolution=resolution)
    Idents(obj_sub) <- detect_cluster_col(obj_sub)
    
    sub_rds <- file.path(out_subdir, sprintf("obj_lineage_%s.rds", gsub("[^A-Za-z0-9]+","_", lin)))
    saveRDS(obj_sub, sub_rds)
    subset_paths[[lin]] <- sub_rds
    
    p_sub <- DimPlot(obj_sub, reduction="umap", group.by=detect_cluster_col(obj_sub),
                     label=TRUE, label.size=UMAP_LABEL_SIZE,
                     pt.size=UMAP_POINT_SIZE, repel=TRUE, raster=TRUE) +
      ggtitle(sprintf("Subset UMAP — %s (new clusters)", lin)) + theme_minimal()
    ggsave(file.path(out_subdir, "UMAP_subset_clusters.png"),
           p_sub, width=PLOT_WIDTH, height=PLOT_HEIGHT, dpi=DPI)
    
    feature_umaps(obj_sub, paste0("subset_", gsub("[^A-Za-z0-9]+","_", lin)), out_subdir,
                  UMAP_POINT_SIZE=UMAP_POINT_SIZE, PLOT_WIDTH=PLOT_WIDTH, PLOT_HEIGHT=PLOT_HEIGHT, DPI=DPI)
    violin_core(obj_sub, paste0("subset_", gsub("[^A-Za-z0-9]+","_", lin)), out_subdir,
                UMAP_LABEL_SIZE=UMAP_LABEL_SIZE, PLOT_WIDTH=PLOT_WIDTH, PLOT_HEIGHT=PLOT_HEIGHT, DPI=DPI)
    
    subset_umaps[[lin]] <- p_sub
  }
  
  section("FOUR-PANEL UMAP")
  choose_lin <- function(pref, available) {
    keep <- intersect(pref, available)
    if (length(keep) >= 2) return(keep[1:2])
    if (length(keep) == 1) {
      others <- setdiff(available, keep)
      return(c(keep, head(others, 1)))
    }
    head(available, 2)
  }
  preferred <- c("Haemopoietic","Non-Haemopoietic")
  available <- names(subset_umaps)
  chosen <- if (length(available) >= 2) choose_lin(preferred, available) else available
  
  p1 <- p_global_clust + ggtitle("Global (original) — Clusters")
  p2 <- p_global_lineage + ggtitle("Global (original) — Lineage")
  p3 <- if (length(chosen) >= 1 && chosen[1] %in% available) subset_umaps[[ chosen[1] ]] + ggtitle(paste0("Subset — ", chosen[1])) else ggplot() + ggtitle("N/A")
  p4 <- if (length(chosen) >= 2 && chosen[2] %in% available) subset_umaps[[ chosen[2] ]] + ggtitle(paste0("Subset — ", chosen[2])) else ggplot() + ggtitle("N/A")
  
  panel_4 <- (p1 | p2) / (p3 | p4)
  ggsave(file.path(out_dir, "UMAP_panel_4.png"), panel_4, width=10*1.8, height=8*1.6, dpi=DPI)
  
  return(invisible(list(
    obj = obj,
    subset_paths = unname(unlist(subset_paths)),
    out_dir = out_dir
  )))
}

## ----------------------------- STEP 6 ----------------------------------------
process_seurat_object_step6 <- function(obj, file_name, output_dir,
                                        dims6 = 30, resolution6 = "0.3,0.5,0.8,1.0",
                                        cluster_doublet_threshold = 70) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  file_base <- gsub("\\.rds$", "", file_name)
  file_base <- gsub("_with_doublets$", "", file_base)
  
  logmsg("==============================================")
  logmsg("Processing (Step 6): %s", file_name)
  logmsg("==============================================")
  
  seurat_v5 <- packageVersion("Seurat") >= "5.0.0"
  if (seurat_v5) logmsg("Detected Seurat v%s", packageVersion("Seurat"))
  
  logmsg("Loaded: %d cells, %d genes", ncol(obj), nrow(obj))
  
  file_output_dir <- file.path(output_dir, file_base)
  dir.create(file_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  obj@meta.data$DoubletDetection <- obj@meta.data$DoubletDetection_Droplet
  
  if (!"DoubletDetection" %in% colnames(obj@meta.data)) {
    logmsg("DoubletDetection column not found in %s — skipping Step 6 for this object.", file_name)
    output_file <- file.path(output_dir, paste0(file_base, "_no_doublet_field.rds"))
    saveRDS(obj, output_file)
    logmsg("Saved (no doublet field): %s\n", output_file)
    return(invisible(NULL))
  }
  
  logmsg("Standardizing DoubletDetection column...")
  if (is.character(obj$DoubletDetection)) {
    logmsg("  Converting character DoubletDetection to logical...")
    obj$DoubletDetection <- tolower(obj$DoubletDetection) == "doublet"
  } else if (is.factor(obj$DoubletDetection)) {
    logmsg("  Converting factor DoubletDetection to logical...")
    obj$DoubletDetection <- tolower(as.character(obj$DoubletDetection)) == "doublet"
  } else if (!is.logical(obj$DoubletDetection)) {
    obj$DoubletDetection <- as.logical(obj$DoubletDetection)
  }
  if (any(is.na(obj$DoubletDetection))) {
    logmsg("  Warning: Found %d NA values in DoubletDetection, treating as singlets",
           sum(is.na(obj$DoubletDetection)))
    obj$DoubletDetection[is.na(obj$DoubletDetection)] <- FALSE
  }
  
  logmsg("Backing up original UMAP and cluster assignments...")
  if ("umap" %in% names(obj@reductions)) {
    obj[["umap_b4_doublet"]] <- obj[["umap"]]
  } else {
    logmsg("Warning: No UMAP found, skipping backup")
  }
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    obj$seurat_clusters_b4_doublet <- obj$seurat_clusters
    obj$idents_b4_doublet <- Idents(obj)
  }
  
  n_doublets <- sum(obj$DoubletDetection, na.rm = TRUE)
  n_total <- ncol(obj)
  pct_doublets <- 100 * n_doublets / n_total
  
  logmsg("DoubletDetection results: %d doublets of %d cells (%.1f%%)",
         n_doublets, n_total, pct_doublets)
  
  if (n_doublets == 0) {
    logmsg("No doublets detected! Saving object as-is.")
    output_file <- file.path(output_dir, paste0(file_base, "_no_doublets.rds"))
    saveRDS(obj, output_file)
    logmsg("Saved: %s\n", output_file)
    return(invisible(NULL))
  }
  
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    cluster_stats <- obj@meta.data %>%
      group_by(seurat_clusters) %>%
      summarise(
        n_cells = n(),
        n_doublets = sum(DoubletDetection, na.rm = TRUE),
        pct_doublets = 100 * n_doublets / n_cells,
        .groups = "drop"
      ) %>% arrange(desc(pct_doublets))
    
    clusters_to_remove <- cluster_stats %>%
      filter(pct_doublets >= cluster_doublet_threshold) %>%
      pull(seurat_clusters)
    
    cluster_stats <- cluster_stats %>%
      mutate(
        action = ifelse(pct_doublets >= cluster_doublet_threshold,
                        paste0("REMOVED_ENTIRE_CLUSTER (>=", cluster_doublet_threshold, "% doublets)"),
                        "Individual doublets removed")
      )
    
    print(cluster_stats)
    write.csv(cluster_stats,
              file.path(file_output_dir, "doublet_rates_by_cluster.csv"),
              row.names = FALSE)
    
    if (length(clusters_to_remove) > 0) {
      logmsg("Clusters with >=%d%% doublet rate (will be completely removed): %s",
             cluster_doublet_threshold, paste(clusters_to_remove, collapse = ", "))
    }
  }
  
  logmsg("Creating pre-removal visualizations...")
  pdf(file.path(file_output_dir, "doublets_before_removal.pdf"),
      width = 18, height = 14)
  
  if ("umap_b4_doublet" %in% names(obj@reductions)) {
    p1 <- DimPlot(obj, reduction = "umap_b4_doublet",
                  group.by = "DoubletDetection",
                  cols = c("FALSE" = "lightgray", "TRUE" = "red")) +
      ggtitle(paste0(file_base, ": DoubletDetection calls")) +
      labs(color = "Doublet")
    
    if ("seurat_clusters_b4_doublet" %in% colnames(obj@meta.data)) {
      p2 <- DimPlot(obj, reduction = "umap_b4_doublet",
                    label = TRUE,
                    group.by = "seurat_clusters_b4_doublet") +
        ggtitle("Original clusters") +
        NoLegend()
      print(p1 | p2)
    } else {
      print(p1)
    }
    
    if (exists("cluster_stats")) {
      p3 <- ggplot(cluster_stats, aes(x = seurat_clusters, y = pct_doublets,
                                      fill = pct_doublets >= cluster_doublet_threshold)) +
        geom_col() +
        geom_hline(yintercept = cluster_doublet_threshold, 
                   linetype = "dashed", color = "red", linewidth = 1) +
        geom_text(aes(label = sprintf("%.1f%%", pct_doublets)), 
                  vjust = -0.5, size = 3) +
        scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "steelblue"),
                          labels = c("TRUE" = "Will remove cluster", 
                                     "FALSE" = "Keep cluster"),
                          name = "") +
        theme_minimal() +
        labs(title = paste0("Doublet rate by cluster (threshold = ", 
                            cluster_doublet_threshold, "%)"),
             x = "Cluster", y = "% Doublets") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      print(p3)
      
      if (length(clusters_to_remove) > 0) {
        obj$cluster_will_remove <- ifelse(obj$seurat_clusters %in% clusters_to_remove,
                                          "Remove", "Keep")
        p4 <- DimPlot(obj, reduction = "umap_b4_doublet",
                      group.by = "cluster_will_remove",
                      cols = c("Remove" = "darkred", "Keep" = "lightgray")) +
          ggtitle(paste0("Clusters to remove (n=", length(clusters_to_remove), ")"))
        print(p4)
      }
    }
  }
  
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    obj$doublet_numeric <- as.numeric(obj$DoubletDetection)
    
    p5 <- ggplot(obj@meta.data, 
                 aes(x = seurat_clusters, y = doublet_numeric, 
                     fill = seurat_clusters)) +
      geom_violin(scale = "width") +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      theme_minimal() +
      labs(title = "Doublet distribution per cluster",
           x = "Cluster", y = "Is Doublet (0=No, 1=Yes)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    print(p5)
    
    cluster_counts <- obj@meta.data %>%
      group_by(seurat_clusters, DoubletDetection) %>%
      summarise(count = n(), .groups = "drop")
    
    p6 <- ggplot(cluster_counts, 
                 aes(x = seurat_clusters, y = count, fill = DoubletDetection)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"),
                        labels = c("FALSE" = "Singlet", "TRUE" = "Doublet")) +
      theme_minimal() +
      labs(title = "Cell counts per cluster",
           x = "Cluster", y = "Number of cells", fill = "") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p6)
  }
  
  dev.off()
  
  logmsg("Removing doublets...")
  cells_to_keep <- rep(TRUE, ncol(obj))
  
  if (exists("clusters_to_remove") && length(clusters_to_remove) > 0) {
    cells_in_bad_clusters <- obj$seurat_clusters %in% clusters_to_remove
    cells_to_keep <- cells_to_keep & !cells_in_bad_clusters
    n_cluster_removed <- sum(cells_in_bad_clusters)
    logmsg("  Removing %d cells from %d high-doublet clusters (>=%d%% doublets)",
           n_cluster_removed, length(clusters_to_remove), cluster_doublet_threshold)
  }
  
  cells_to_keep <- cells_to_keep & !obj$DoubletDetection
  obj_clean <- obj[, cells_to_keep]
  
  n_removed <- ncol(obj) - ncol(obj_clean)
  logmsg("After doublet removal: %d cells", ncol(obj_clean))
  logmsg("Removed: %d cells (%.1f%%)", n_removed, 100 * n_removed / ncol(obj))
  
  logmsg("Re-clustering after doublet removal...")
  
  if ("RNA" %in% names(obj_clean@assays)) {
    if (inherits(obj_clean[["RNA"]], "Assay5")) {
      logmsg("Detected Seurat v5 - joining layers...")
      obj_clean[["RNA"]] <- JoinLayers(obj_clean[["RNA"]])
    }
  }
  
  if ("RNA" %in% names(obj_clean@assays)) {
    if ("scale.data" %in% slotNames(obj_clean[["RNA"]])) {
      slot(obj_clean[["RNA"]], "scale.data") <- matrix()
    }
  }
  
  obj_clean <- FindVariableFeatures(obj_clean,
                                    selection.method = "vst",
                                    nfeatures = 2000,
                                    verbose = FALSE)
  
  obj_clean <- ScaleData(obj_clean, verbose = FALSE)
  obj_clean <- RunPCA(obj_clean, npcs = 50, verbose = FALSE)
  obj_clean <- FindNeighbors(obj_clean, dims = 1:dims6, verbose = FALSE)
  
  resolutions <- as.numeric(unlist(strsplit(resolution6, ",")))
  obj_clean <- FindClusters(obj_clean, resolution = resolutions, verbose = FALSE)
  
  default_res <- resolutions[ceiling(length(resolutions)/2)]
  res_col <- paste0("RNA_snn_res.", default_res)
  if (res_col %in% colnames(obj_clean@meta.data)) {
    Idents(obj_clean) <- obj_clean@meta.data[[res_col]]
    obj_clean$seurat_clusters <- as.factor(obj_clean@meta.data[[res_col]])
  }
  
  obj_clean <- RunUMAP(obj_clean, dims = 1:dims6, verbose = FALSE)
  
  logmsg("New clustering complete! New clusters: %d (at resolution %.1f)",
         length(unique(obj_clean$seurat_clusters)), default_res)
  
  logmsg("Creating before/after comparison plots...")
  pdf(file.path(file_output_dir, "before_after_doublet_removal.pdf"),
      width = 20, height = 16)
  
  if ("umap_b4_doublet" %in% names(obj@reductions) &&
      "seurat_clusters_b4_doublet" %in% colnames(obj@meta.data)) {
    
    p_before_1 <- DimPlot(obj, reduction = "umap_b4_doublet",
                          group.by = "seurat_clusters_b4_doublet",
                          label = TRUE, label.size = 5) +
      ggtitle(sprintf("BEFORE: %d cells, %d clusters",
                      ncol(obj),
                      length(unique(obj$seurat_clusters_b4_doublet)))) +
      NoLegend()
    
    p_before_2 <- DimPlot(obj, reduction = "umap_b4_doublet",
                          group.by = "DoubletDetection",
                          cols = c("FALSE" = "lightgray", "TRUE" = "red")) +
      ggtitle(sprintf("Doublets marked (%d cells, %.1f%%)",
                      n_doublets, pct_doublets))
    
    print(p_before_1 | p_before_2)
  }
  
  p_after_1 <- DimPlot(obj_clean, reduction = "umap",
                       label = TRUE, label.size = 5) +
    ggtitle(sprintf("AFTER: %d cells, %d clusters",
                    ncol(obj_clean),
                    length(unique(obj_clean$seurat_clusters)))) +
    NoLegend()
  
  p_after_2 <- DimPlot(obj_clean, reduction = "umap") +
    ggtitle("Cleaned data (all doublets removed)")
  
  print(p_after_1 | p_after_2)
  
  dev.off()
  
  output_file <- file.path(output_dir, paste0(file_base, "_cleaned.rds"))
  saveRDS(obj_clean, output_file)
  logmsg("Saved: %s", output_file)
  invisible(obj_clean)
}

run_step6_on_dir <- function(input_dir, outdir_step6, dims6, resolution6, cluster_doublet_threshold) {
  rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
  rds_files <- rds_files[grepl("obj_lineage_", basename(rds_files))]  # limit to per-lineage outputs
  
  if (length(rds_files) == 0) {
    logmsg("No per-lineage RDS files found under: %s", input_dir)
    return(invisible(NULL))
  }
  
  logmsg("Found %d RDS file(s) to process in Step 6:", length(rds_files))
  for (f in rds_files) logmsg("  - %s", basename(f))
  
  dir.create(outdir_step6, showWarnings = FALSE, recursive = TRUE)
  
  for (rds in rds_files) {
    obj <- readRDS(rds)
    if (!inherits(obj, "Seurat")) {
      logmsg("Skipping %s - not a Seurat object", basename(rds))
      next
    }
    tryCatch({
      process_seurat_object_step6(obj,
                                  file_name = basename(rds),
                                  output_dir = outdir_step6,
                                  dims6 = dims6,
                                  resolution6 = resolution6,
                                  cluster_doublet_threshold = cluster_doublet_threshold)
    }, error = function(e) {
      logmsg("ERROR processing %s: %s", basename(rds), e$message)
      logmsg("Traceback follows:")
      print(traceback())
    })
  }
  logmsg("Step 6 completed for directory: %s", normalizePath(outdir_step6))
}

## ------------------------------ MAIN -----------------------------------------
main <- function() {
  ## Parse umap dims string to a sequence (e.g., "1:30" -> 1:30)
  umap_dims <- tryCatch({
    eval(parse(text = opt$umap_dims))
  }, error = function(e) {
    warning("Invalid --umap_dims '", opt$umap_dims, "'. Falling back to 1:30")
    1:30
  })
  
  step5_res <- run_step5(
    in_rds     = opt$in_rds,
    annot_csv  = opt$annot_csv,
    out_dir    = opt$outdir_step5,
    assay      = opt$assay,
    use_sct    = opt$use_sct,
    n_pcs      = opt$n_pcs,
    umap_dims  = umap_dims,
    resolution = opt$resolution5
  )
  
  if (isTRUE(opt$run_step6)) {
    section("RUN STEP 6 ON STEP5 PER-LINEAGE OBJECTS")
    step5_subsets_dir <- file.path(opt$outdir_step5, "subsets")
    run_step6_on_dir(
      input_dir = step5_subsets_dir,
      outdir_step6 = opt$outdir_step6,
      dims6 = opt$dims6,
      resolution6 = opt$resolution6,
      cluster_doublet_threshold = opt$cluster_doublet_threshold
    )
  } else {
    logmsg("Skipping Step 6 (doublet removal) as requested.")
  }
  
  section("DONE")
}

if (identical(environment(), globalenv())) {
  main()
}



#Rscript v2integrated_lineage_doublet_split_panGI_pipeline.R \
#--in_rds ./MERGED_PIPELINE_OUTPUT_V4/object_final_annotated_panGI.rds \
#--outdir_step5  LINEAGE_DOUBLET_ANALYSIS_PANGI_V7 \
#--run_step6 TRUE \
#--split_by pangi


#Rscript v2integrated_lineage_doublet_split_panGI_pipeline.R \
#--in_rds ./MERGED_PIPELINE_OUTPUT_V4/object_final_annotated_panGI.rds \
#--outdir_step5  LINEAGE_DOUBLET_ANALYSIS_PANGI_V8 \
#--run_step6 TRUE \
#--split_by csv