#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mem=13Gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=glmeeks@ucdavis.edu
#SBATCH --array=1-100   # Number of iterations in the loop
#SBATCH --partition=production
#SBATCH --job-name=fake_meqtl_meta
#SBATCH --output=output_%A_%a.txt
#SBATCH --error=error_%A_%a.txt

source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate meta_R

echo "$SLURM_ARRAY_TASK_ID"
pop="$1"
# Check if Conda environment was activated successfully
if [ $? -eq 0 ]; then
    echo "Conda environment activated successfully."
else
    echo "Error: Failed to activate Conda environment ."
    exit 1
fi


# Run R script for each iteration
echo "Rscript metagen_metanalyis.R $SLURM_ARRAY_TASK_ID"



Rscript metagen_metanalysis.R $SLURM_ARRAY_TASK_ID 

