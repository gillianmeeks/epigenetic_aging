#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=25:00:00
#SBATCH --mem=35Gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=glmeeks@ucdavis.edu
#SBATCH --array=1-100   # Number of iterations in the loop
#SBATCH --partition=production
#SBATCH --job-name=notheritable
#SBATCH --output=output_%A_%a.txt
#SBATCH --error=error_%A_%a.txt

source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh

conda activate R_4.1.3

echo "$SLURM_ARRAY_TASK_ID"
task="$SLURM_ARRAY_TASK_ID"
Rscript LOOCV_train_test_split_cv.R "$task"


