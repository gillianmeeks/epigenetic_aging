#!/bin/bash/

results_dir="/share/hennlab/users/glmeeks/age_methylation/FUSION/FUSION_results_FINAL"
data_dir="/share/hennlab/users/glmeeks/age_methylation/FUSION/FUSION_data/FUSION_pheno"

source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate FUSION_env

for P in KHS Himba Baka;do
  # Create position file
  python gen_pos_all_hsq.py ${data_dir}/${P}/${P}_hg19_meth.bed.gz ${results_dir}/${P} ${P} > ${results_dir}/${P}_hsq_pos_file4.tsv
  # Extract R2, hsq
  #Rscript FUSION.assoc_test_gutted_SSG.R --out ${results_dir}/${P}_R2_hsq_all_4.tsv --weights ${results_dir}/${P}_hsq_pos_file4.tsv --weights_dir ${results_dir}/${P}
done
#test
