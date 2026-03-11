library(Seurat)
library(dplyr)

obj <- readRDS("/gpfs/gsfs12/users/cheungf/ACEPD/LINEAGE_DOUBLET_ANALYSIS_PANGI_V9/object_final_annotated_panGI_withDay.rds")

run_markers <- function(
    obj,
    ident_col,
    out_prefix,
    only.pos = TRUE,
    test.use = "wilcox",
    min.pct = 0.25,
    logfc.threshold = 0.25,
    padj_cutoff = 0.05,
    pct1_cutoff = 0.3,
    top_n = 10,
    min_cells_per_ident = NULL
) {
  # sanity check
  if (!ident_col %in% colnames(obj@meta.data)) {
    stop(paste0("Identity column not found in meta.data: ", ident_col))
  }
  
  Idents(obj) <- ident_col
  
  # Optional: drop small identity groups (useful for seurat_clusters small clusters)
  if (!is.null(min_cells_per_ident)) {
    id_sizes <- table(Idents(obj))
    keep <- names(which(id_sizes >= min_cells_per_ident))
    obj <- subset(obj, idents = keep)
    message(out_prefix, ": kept ", length(keep), " identities with >= ", min_cells_per_ident, " cells")
  }
  
  # Find markers
  m_all <- FindAllMarkers(
    obj,
    only.pos = only.pos,
    test.use = test.use,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold
  )
  
  # Save full table
  write.csv(m_all, paste0(out_prefix, "_ALL.csv"), row.names = FALSE)
  
  # Top N per identity
  m_top <- m_all %>%
    filter(p_val_adj < padj_cutoff, pct.1 > pct1_cutoff) %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = top_n, with_ties = FALSE) %>%
    ungroup()
  
  write.csv(m_top, paste0(out_prefix, "_top", top_n, ".csv"), row.names = FALSE)
  
  return(list(all = m_all, top = m_top))
}

# 2) Markers by final UMAP clusters (seurat_clusters)
res_seurat_clusters <- run_markers(
  obj = obj,
  ident_col = "seurat_clusters",
  out_prefix = "V1_markers_seurat_clusters",
  min_cells_per_ident = 200     # optional; adjust or set to NULL to keep all clusters
)


# 1) Markers by Pan-GI L2 annotation
res_panGI_L2 <- run_markers(
  obj = obj,
  ident_col = "PanGI_L2",
  out_prefix = "v1_markers_PanGI_L2",
  min_cells_per_ident = NULL   # usually keep all cell types
)


message("Done. Wrote CSVs for PanGI_L2 and seurat_clusters marker sets.")
