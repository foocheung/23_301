#!/bin/bash
#SBATCH -J Step1DV3
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH -t 48:00:00

set -euo pipefail
mkdir -p logs

### Load your site’s R + pandoc (names may differ on your cluster)
#module purge

### Knit the report (HTML). Change filename to your Rmd.
quarto render Step1D_c.qmd 
