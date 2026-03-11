#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(escape)
  library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)

# Usage:
# Rscript escape_chunk.R <seurat_rds> <outdir> <chunk_size> <task_id>
seurat_rds <- args[1]
outdir     <- args[2]
chunk_size <- as.integer(args[3])
task_id    <- as.integer(args[4])

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("Reading Seurat: ", seurat_rds)
seu <- readRDS(seurat_rds)
stopifnot("PanGI_L2" %in% colnames(seu@meta.data))

DefaultAssay(seu) <- "RNA"

cells_all <- colnames(seu)
n_cells   <- length(cells_all)

# Define chunk indices
n_chunks <- ceiling(n_cells / chunk_size)
stopifnot(task_id >= 1, task_id <= n_chunks)

start_idx <- (task_id - 1) * chunk_size + 1
end_idx   <- min(task_id * chunk_size, n_cells)

cells_i <- cells_all[start_idx:end_idx]
message("Task ", task_id, "/", n_chunks, " cells [", start_idx, ":", end_idx, "] n=", length(cells_i))

# ---------------------------
# GO:BP gene sets (msigdbr >= 10)
# ---------------------------
go_df <- msigdbr(
  species       = "Homo sapiens",
  collection    = "C5",
  subcollection = "BP"
)[, c("gs_name", "gene_symbol")]
go_df <- unique(go_df)

go_sets <- split(go_df$gene_symbol, go_df$gs_name)

# Expression for this chunk (genes x cells)
expr_i <- GetAssayData(seu, assay = "RNA", layer = "data")[, cells_i, drop = FALSE]
genes_present <- rownames(expr_i)

# ---------------------------
# IMPORTANT: filter huge pathways so UCell won't crash
# ---------------------------
# UCell default maxRank is typically 1500; keep signatures <= 1500 to be safe
MIN_SET <- 10
MAX_SET <- 1500

go_sets_filt <- lapply(go_sets, function(g) intersect(g, genes_present))
go_sets_filt <- go_sets_filt[lengths(go_sets_filt) >= MIN_SET & lengths(go_sets_filt) <= MAX_SET]

message("GO pathways loaded: ", length(go_sets),
        " | after filtering: ", length(go_sets_filt),
        " (kept size ", MIN_SET, "–", MAX_SET, " after intersect)")

if (length(go_sets_filt) == 0) {
  stop("No GO:BP pathways left after filtering. Try increasing MAX_SET or check gene symbols.")
}

# ---------------------------
# Run escape (UCell)
# ---------------------------
sc_i <- escape.matrix(
  input.data = expr_i,
  gene.sets  = go_sets_filt,
  method     = "UCell"
)

# Ensure pathways x cells and correct cell order
# (escape outputs can vary; this forces a consistent pathways-by-cells matrix)
if (is.null(colnames(sc_i)) || is.null(rownames(sc_i))) {
  stop("escape returned matrix missing dimnames.")
}
if (!identical(colnames(sc_i), cells_i) && identical(rownames(sc_i), cells_i)) {
  sc_i <- t(sc_i)
}
sc_i <- sc_i[, cells_i, drop = FALSE]

# Save chunk scores + the cell list for safety
chunk_file <- file.path(outdir, sprintf("escape_scores_chunk_%05d.rds", task_id))
meta_file  <- file.path(outdir, sprintf("cells_chunk_%05d.rds", task_id))

saveRDS(as.matrix(sc_i), chunk_file)
saveRDS(cells_i, meta_file)

message("Wrote: ", chunk_file)
message("Wrote: ", meta_file)
