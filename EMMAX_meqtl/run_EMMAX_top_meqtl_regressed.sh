#just feed in the regressed beta_matrix instead

#creat BN matrix
#plink --bfile [bed_prefix] (or --file [ped_prefix]) --recode12 --output-missing-genotype 0 --transpose --out [tped_prefix]
#/share/hennlab/progs/emmax-beta-07Mar2010/emmax-kin -v -d 10 all_pops_merged_segregating_pos_ids_maf_01

###MEQTL regressed ewas####


/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t testbeta_meqtl_regress_KHS -p age_KHS.txt -k KHS_95_MAC02_pos_id.BN.kinf -c covs_KHS_sex_batch_PC13_meth_Neutro_genoPCs.txt -Z -o KHS_sex_batch_PC13_meth_Neutro_genoPCs_meqtl_regressed



/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t testbeta_meqtl_regress_Himba -p age_Himba.txt -k Himba_95_MAC02_pos_id.BN.kinf -c  covs_Himba_sex_cpPCs_Neutro_genoPCs.txt -Z -o Himba_sex_cpPCs_Neutro_genoPCs_meqtl_regressed




/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t testbeta_meqtl_regress_Baka -p age_Baka.txt -k Baka_95_MAC02_pos_id.BN.kinf -c covs_Baka_sex_PC23_meth_Neutro_genoPC12345.txt -Z -o Baka_sex_methPC23_Neutro_genoPC12345_meqtl_regressed


/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t testbeta_meqtl_regress_all -p age_all.txt -k all_pops_merged_segregating_pos_ids_maf_01.BN.kinf -c covs_all.txt -Z -o all_covs_meqtl_regressed


############################STANDARD EWAS###

/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_no_regress_KHS -p age_KHS.txt -k KHS_95_MAC02_pos_id.BN.kinf -c covs_KHS_sex_batch_PC13_meth_Neutro_genoPCs.txt -Z -o KHS_sex_batch_PC13_meth_Neutro_genoPCs_no_meqtl_regress



/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_no_regress_Himba -p age_Himba.txt -k Himba_95_MAC02_pos_id.BN.kinf -c  covs_Himba_sex_cpPCs_Neutro_genoPCs.txt -Z -o Himba_sex_cpPCs_Neutro_genoPCs_no_meqtl_regress



/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_meqtl_regress_Baka -p age_Baka.txt -k Baka_95_MAC02_pos_id.BN.kinf -c covs_Baka_sex_PC23_meth_Neutro_genoPC12345.txt -Z -o Baka_sex_methPC23_Neutro_genoPC12345_no_meqtl_regress


/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_all_no_regress -p age_all.txt -k all_pops_merged_segregating_pos_ids_maf_01.BN.kinf -c covs_all.txt -Z -o all_covs_no_meqtl_regressed


#################FUSION REGRESSED######
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t test_beta_meqtl_regress_fusionall -p age_all.txt -k all_pops_merged_segregating_pos_ids_maf_01.BN.kinf -c covs_all.txt -Z -o all_covs_fusion_regressed


/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t test_beta_meqtl_regress_fusionKHS -p age_KHS.txt -k KHS_95_MAC02_pos_id.BN.kinf -c covs_KHS_sex_batch_PC13_meth_Neutro_genoPCs.txt -Z -o KHS_sex_batch_PC13_meth_Neutro_genoPCs_fusion_regressed


/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t test_beta_meqtl_regress_fusionBaka -p age_Baka.txt -k  Baka_95_MAC02_pos_id.BN.kinf -c covs_Baka_sex_PC23_meth_Neutro_genoPC12345.txt -Z -o Baka_sex_methPC23_Neutro_genoPC12345_fusion_regressed

/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t test_beta_meqtl_regress_fusionHimba -p age_Himba.txt -k Himba_95_MAC02_pos_id.BN.kinf -c  covs_Himba_sex_cpPCs_Neutro_genoPCs.txt -Z -o Himba_sex_cpPCs_Neutro_genoPCs_fusion_regressed
