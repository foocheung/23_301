#!/bin/bash
#SBATCH --job-name=escape_ucell
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --array=1-29
#SBATCH --output=logs/escape_%A_%a.out
#SBATCH --error=logs/escape_%A_%a.err

#module load R  # adjust for your cluster

SEURAT_RDS="object_final_annotated_panGI_withDay.rds"
OUTDIR="escape_chunks"
CHUNK_SIZE=5000

mkdir -p logs
mkdir -p "${OUTDIR}"

Rscript escape_chunk_GO.R "${SEURAT_RDS}" "${OUTDIR}" "${CHUNK_SIZE}" "${SLURM_ARRAY_TASK_ID}"
##Rscript -e 'seu<-readRDS("object_final_annotated_panGI_withDay.rds"); cat(ceiling(ncol(seu)/5000), "\n")'
##Then replace --array=1-<NCHUNKS> with that number.