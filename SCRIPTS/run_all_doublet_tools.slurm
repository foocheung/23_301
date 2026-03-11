#!/bin/bash
#SBATCH --job-name=ALL_tools
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --output=logs/W_%x_%j.out
#SBATCH --error=logs/W_%x_%j.err


module load r


# Run your hard-coded R script (no args)
bash run_all_doublet_tools.sh 


