#!/bin/bash

# Define function to execute EMMAX command for a specific iteration
run_emmax() {
    i="$1"
    command="$2"
    echo "Running EMMAX for iteration $i"
    $command
}

# Define the EMMAX commands for each population
commands=(


  "/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t fake_perms/{i}beta_fake_meqtl_regress_KHS -p age_KHS.txt -k KHS_95_MAC02_pos_id.BN.kinf -c covs_KHS_sex_batch_PC13_meth_Neutro_genoPCs.txt -Z -o fake_perms/{i}beta_fake_meqtl_regress_KHS"
  
  "/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t fake_perms/{i}beta_fake_meqtl_regress_Himba -p age_Himba.txt -k Himba_95_MAC02_pos_id.BN.kinf -c covs_Himba_sex_cpPCs_Neutro_genoPCs.txt -Z -o fake_perms/{i}beta_fake_meqtl_regress_Himba"
  
  "/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t fake_perms/{i}beta_fake_meqtl_regress_Baka -p age_Baka.txt -k Baka_95_MAC02_pos_id.BN.kinf -c covs_Baka_sex_PC23_meth_Neutro_genoPC12345.txt  -Z -o fake_perms/{i}beta_fake_meqtl_regress_Baka"
  


)

# Iterate through iterations 1 to 100
for i in {1..100}; do
    # Run each EMMAX command in parallel
    for command in "${commands[@]}"; do
        # Run EMMAX command for current iteration in background
        run_emmax "$i" "${command//\{i\}/$i}" &
        # Limit parallelism to 30 cores
        if (( $(jobs -r -p | wc -l) >= 30 )); then
            wait -n
        fi
    done
done

# Wait for all background jobs to finish
wait




# 
# 


# 


#  "/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t fake_perms/{i}beta_fake_meqtl_regress_all -p age_all.txt -k all_pops_merged_segregating_pos_ids_maf_01.BN.kinf -c covs_all.txt -Z -o fake_perms/{i}all_covs_meqtl_regressed"
