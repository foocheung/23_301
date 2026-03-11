#!/bin/bash
#SBATCH -J demux_batch_filt
#SBATCH -o logs/demux_batch_filt.%j.out
#SBATCH -e logs/demux_batch_filt.%j.err
#SBATCH -t 48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G

set -euo pipefail
 

#module purge

module load r

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

Rscript /data/chi/PROJECTS/23-301_Lionakis_APECED/Bioinformatics/Results/22NWNCLT3/Flowcell_22NWNCLT3/cellranger_outs/run_dropulation.R

