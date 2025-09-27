#!/bin/bash
#SBATCH --job-name=divalike_bears
#SBATCH --account=bisc033844
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=divalike_run_%j.out
#SBATCH --error=divalike_run_%j.err

# Activate conda env via Miniforge
source /software/local/languages/miniforge3/etc/profile.d/conda.sh
conda activate bears

cd /user/work/ay21702/mitogenome_project/biogeo/run2

# Run the R batch script
Rscript divalike_batch.R

