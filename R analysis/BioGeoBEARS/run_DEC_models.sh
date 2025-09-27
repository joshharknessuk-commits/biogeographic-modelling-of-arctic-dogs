#!/bin/bash
#SBATCH --job-name=DEC_models
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --account=bisc033844
#SBATCH --output=DEC_models.%j.out
#SBATCH --error=DEC_models.%j.err

# Load conda from Miniforge
source /software/local/languages/miniforge3/etc/profile.d/conda.sh
conda activate bears

cd "$SLURM_SUBMIT_DIR"

Rscript dec_batch.R

