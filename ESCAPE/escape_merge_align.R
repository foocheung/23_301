#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)

# Usage:
# Rscript escape_merge_align.R <seurat_rds> <chunk_outdir> <out_seurat_rds> <out_scores_rds> <n_chunks> <fill_value>
# fill_value: "NA" or "0"
seurat_rds     <- args[1]
chunk_outdir   <- args[2]
out_seurat_rds <- args[3]
out_scores_rds <- args[4]
n_chunks       <- as.integer(args[5])
fill_value     <- if (length(args) >= 6) args[6] else "NA"

fill_is_na <- identical(toupper(fill_value), "NA")

message("Reading Seurat: ", seurat_rds)
seu <- readRDS(seurat_rds)
cells_all <- colnames(seu)

# ---------------------------
# 1) Load all chunks
# ---------------------------
chunks <- vector("list", n_chunks)
for (i in seq_len(n_chunks)) {
  f <- file.path(chunk_outdir, sprintf("escape_scores_chunk_%05d.rds", i))
  if (!file.exists(f)) stop("Missing chunk file: ", f)
  m <- readRDS(f)
  chunks[[i]] <- as.matrix(m)
  message("Loaded chunk ", i, " dim=", paste(dim(chunks[[i]]), collapse="x"))
}

# ---------------------------
# 2) Create union of pathways (row names)
# ---------------------------
all_pathways <- unique(unlist(lapply(chunks, rownames)))
message("Union pathways: ", length(all_pathways))

# ---------------------------
# 3) Reindex each chunk to union pathways (fill missing)
# ---------------------------
reindexed <- vector("list", n_chunks)

for (i in seq_len(n_chunks)) {
  m <- chunks[[i]]
  # Ensure colnames exist
  if (is.null(colnames(m))) stop("Chunk ", i, " has NULL colnames.")
  # Create template matrix
  out <- matrix(
    if (fill_is_na) NA_real_ else 0,
    nrow = length(all_pathways),
    ncol = ncol(m),
    dimnames = list(all_pathways, colnames(m))
  )
  # Fill overlapping pathways
  idx <- match(rownames(m), all_pathways)
  out[idx, ] <- m
  reindexed[[i]] <- out
  rm(out, m); gc()
}

# ---------------------------
# 4) Combine chunks
# ---------------------------
escape_scores <- do.call(cbind, reindexed)

# Force final column order to match Seurat exactly
common <- intersect(colnames(escape_scores), cells_all)
if (length(common) != length(cells_all)) {
  stop("Merged matrix does not contain all Seurat cells. Have=", length(common), " Need=", length(cells_all))
}
escape_scores <- escape_scores[, cells_all, drop=FALSE]

message("Final merged dims: ", paste(dim(escape_scores), collapse="x"))

# ---------------------------
# 5) Add assay + save
# ---------------------------
seu[["ESCAPE"]] <- CreateAssayObject(data = escape_scores)

saveRDS(seu, out_seurat_rds)
saveRDS(escape_scores, out_scores_rds)

message("Saved Seurat: ", out_seurat_rds)
message("Saved scores: ", out_scores_rds)
