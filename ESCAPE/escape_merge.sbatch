#!/bin/bash
#SBATCH --job-name=escape_merge
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=logs/escape_merge_%j.out
#SBATCH --error=logs/escape_merge_%j.err

module load R  # adjust

SEURAT_RDS="object_final_annotated_panGI_withDay.rds"
CHUNK_OUTDIR="escape_chunks"
OUT_SEURAT="v2_object_final_annotated_panGI_withDay_ESCAPE_Reactome_ALLCELLS_ARRAY.rds"
OUT_SCORES="v2_ESCAPE_Reactome_UCell_scores_ALLCELLS_ARRAY.rds"
NCHUNKS=29

Rscript escape_merge_align.R "${SEURAT_RDS}" "${CHUNK_OUTDIR}" "${OUT_SEURAT}" "${OUT_SCORES}" "${NCHUNKS}"
