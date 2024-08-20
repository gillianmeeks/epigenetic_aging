#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=23:00:00
#SBATCH --mem=12Gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=glmeeks@ucdavis.edu
#SBATCH --array=1-100   # Number of iterations in the loop
#SBATCH --partition=production
#SBATCH --job-name=all_fake_meqtl_EMMAX
#SBATCH --output=output_%A_%a.txt
#SBATCH --error=error_%A_%a.txt


source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh

conda activate meta_R

echo "$SLURM_ARRAY_TASK_ID"

# Check if Conda environment was activated successfully
if [ $? -eq 0 ]; then
    echo "Conda environment activated successfully."
else
    echo "Error: Failed to activate Conda environment ."
    exit 1
fi


cd /share/hennlab/users/glmeeks/age_methylation/age_prediction/

# Run R script for each iteration
echo "Rscript fake_permute_meqtl_emmax_prep.R "$SLURM_ARRAY_TASK_ID"" Himba



Rscript fake_permute_meqtl_emmax_prep.R "$SLURM_ARRAY_TASK_ID" Himba
