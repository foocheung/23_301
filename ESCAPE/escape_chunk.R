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

# Reactome gene sets
reactome_df <- msigdbr(
  species     = "Homo sapiens",
  category    = "C2",
  subcategory = "CP:REACTOME"
)[, c("gs_name", "gene_symbol")]

reactome_df   <- unique(reactome_df)
reactome_sets <- split(reactome_df$gene_symbol, reactome_df$gs_name)
message("Reactome pathways: ", length(reactome_sets))

# Expression for this chunk (genes x cells)
expr_i <- GetAssayData(seu, assay="RNA", layer="data")[, cells_i, drop=FALSE]

# Run escape
sc_i <- escape.matrix(
  input.data = expr_i,
  gene.sets  = reactome_sets,
  method     = "UCell"
)

# Ensure pathways x cells and correct cell order
if (!is.null(colnames(sc_i)) && !identical(colnames(sc_i), cells_i)) sc_i <- t(sc_i)
if (is.null(colnames(sc_i))) stop("escape returned NULL colnames")
sc_i <- sc_i[, cells_i, drop=FALSE]

# Save chunk scores + the cell list for safety
chunk_file <- file.path(outdir, sprintf("escape_scores_chunk_%05d.rds", task_id))
meta_file  <- file.path(outdir, sprintf("cells_chunk_%05d.rds", task_id))

saveRDS(as.matrix(sc_i), chunk_file)
saveRDS(cells_i, meta_file)

message("Wrote: ", chunk_file)
message("Wrote: ", meta_file)
