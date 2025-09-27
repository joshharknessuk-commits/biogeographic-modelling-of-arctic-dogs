#!/bin/bash
#SBATCH --job-name=prefetch
#SBATCH --output=prefetch_%j.out
#SBATCH --error=prefetch_%j.err
#SBATCH --account=bisc033844
#SBATCH --partition=teach_cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00

module load sra-tools

prefetch --max-size 100G --option-file /user/work/ay21702/mitogenome_project/scripts+ref/sra_list.txt --output-directory /user/work/ay21702/mitogenome_project/sra_raw/

