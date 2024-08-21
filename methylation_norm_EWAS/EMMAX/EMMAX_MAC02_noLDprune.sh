
#transpose plink file
plink --bfile [bed] --recode12 --output-missing-genotype 0 --transpose --out [out_prefix]

#make grms
/share/hennlab/progs/emmax-beta-07Mar2010/emmax-kin -v -d 10 [tped]



##do the emmax
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t [beta_tped] -p [age] -k  [kinf]  -c [covar] -Z -o

##-Z is the dosage option


############BAKA#######

/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf  -Z -o Baka_no_covars
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf  -c covs_Baka_sex.txt -Z -o Baka_sex
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf  -c covs_Baka_sex_Neutro.txt -Z -o Baka_sex_Neutro
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf  -c covs_Baka_sex_PC1_meth.txt -Z -o Baka_sex_PC1_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC12_meth.txt -Z -o Baka_sex_PC12_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC13_meth.txt -Z -o Baka_sex_PC13_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC14_meth.txt -Z -o Baka_sex_PC14_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC23_meth.txt -Z -o Baka_sex_PC23_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC24_meth.txt -Z -o Baka_sex_PC24_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC34_meth.txt -Z -o Baka_sex_PC34_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf  -c covs_Baka_sex_PC123_meth.txt -Z -o Baka_sex_PC123_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC1234_meth.txt -Z -o Baka_sex_PC1234_meth


/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf  -c covs_Baka_sex_PC24_meth.txt -Z -o Baka_sex_PC24_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC24_meth_Neutro.txt -Z -o Baka_sex_PC24_meth_Neutro
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf  -c covs_Baka_sex_PC24_meth_Bantu.txt -Z -o Baka_sex_PC24_meth_Bantu
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC24_meth_genoPC12345.txt -Z -o Baka_sex_PC24_meth_genoPC12345
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_ctrl.txt -Z -o Baka_sex_ctrl
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf  -c covs_Baka_sex_ctrl_Neutro.txt -Z -o Baka_sex_ctrl_Neutro


/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC24_meth_Neutro_Bantu_ctrl.txt -Z -o Baka_sex_PC24_meth_Neutro_Bantu_ctrl
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC24_meth_Neutro_genoPC12345_Bantu.txt -Z -o Baka_sex_PC24_meth_Neutro_genoPC12345_Bantu
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf  -c covs_Baka_sex_PC24_meth_Neutro_Bantu.txt -Z -o Baka_sex_PC24_meth_Neutro_Bantu
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC24_meth_Neutro_genoPC12345.txt -Z -o Baka_sex_PC24_meth_Neutro_genoPC12345
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC24_meth_Neutro_genoPC1234.txt -Z -o Baka_sex_PC24_meth_Neutro_genoPC1234
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC24_meth_Neutro_genoPC123.txt -Z -o Baka_sex_PC24_meth_Neutro_genoPC123
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC24_meth_Neutro_genoPC12.txt -Z -o Baka_sex_PC24_meth_Neutro_genoPC12
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC24_meth_Neutro_genoPC1.txt -Z -o Baka_sex_PC24_meth_Neutro_genoPC1

/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC23_meth_Neutro_Bantu_ctrl.txt -Z -o Baka_sex_PC23_meth_Neutro_Bantu_ctrl
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC23_meth_Neutro_Bantu.txt -Z -o Baka_sex_PC23_meth_Neutro_Bantu
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Baka -p age_Baka.txt -k Baka_95_MAC02.BN.kinf -c covs_Baka_sex_PC23_meth_Neutro_genoPC12345.txt -Z -o Baka_sex_PC23_meth_Neutro_genoPC12345


############HIMBA#######

/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Himba -p age_Himba.txt -k Himba_95_MAC02.BN.kinf -Z -o Himba_nocovs
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Himba -p age_Himba.txt -k Himba_95_MAC02.BN.kinf -c covs_Himba_sex_PC1_meth.txt -Z -o Himba_sex_PC1_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Himba -p age_Himba.txt -k Himba_95_MAC02.BN.kinf -c covs_Himba_sex_PC12345_meth.txt -Z -o Himba_sex_PC12345_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Himba -p age_Himba.txt -k Himba_95_MAC02.BN.kinf -c covs_Himba_sex_cpPCs.txt -Z -o Himba_sex_cpPCs
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Himba -p age_Himba.txt -k Himba_95_MAC02.BN.kinf -c covs_Himba_sex_cpPCs_Neutro.txt -Z -o Himba_sex_cpPCs_Neutro
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Himba -p age_Himba.txt -k Himba_95_MAC02.BN.kinf -c covs_Himba_sex_cpPCs_Neutro_genoPCs.txt -Z -o Himba_sex_cpPCs_Neutro_genoPCs
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_Himba -p age_Himba.txt -k Himba_95_MAC02.BN.kinf -c covs_Himba_sex_batch_Neutro_genoPCs.txt -Z -o Himba_sex_batch_Neutro_genoPCs


###########KHS##############

/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -Z -o KHS_nocovs
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC1_meth.txt -Z -o KHS_sex_PC1_meth

/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC12_meth.txt -Z -o KHS_sex_PC12_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC13_meth.txt -Z -o KHS_sex_PC13_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC14_meth.txt -Z -o KHS_sex_PC14_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC23_meth.txt -Z -o KHS_sex_PC23_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC24_meth.txt -Z -o KHS_sex_PC24_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC34_meth.txt -Z -o KHS_sex_PC34_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC123_meth.txt -Z -o KHS_sex_PC123_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC1234_meth.txt -Z -o KHS_sex_PC1234_meth



/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC12345_meth.txt -Z -o KHS_sex_PC12345_meth
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_cpPCs.txt -Z -o KHS_sex_cpPCs
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_cpPCs_Neutro.txt -Z -o KHS_sex_cpPCs_Neutro
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_cpPCs_Neutro_genoPCs.txt -Z -o KHS_sex_cpPCs_Neutro_genoPCs
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_batch_Neutro_genoPCs.txt -Z -o KHS_sex_batch_Neutro_genoPCs
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC1_meth_Neutro_genoPCs.txt -Z -o KHS_sex_PC1_meth_Neutro_genoPCs
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC12345_meth_Neutro_genoPCs.txt -Z -o KHS_sex_PC12345_meth_Neutro_genoPCs



/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC1_meth_cpPCs_Neutro_genoPCs.txt -Z -o KHS_sex_PC1_meth_cpPCs_Neutro_genoPCs
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC12_meth_cpPCs_Neutro_genoPCs.txt -Z -o KHS_sex_PC12_meth_cpPCs_Neutro_genoPCs
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC123_meth_cpPCs_Neutro_genoPCs.txt -Z -o KHS_sex_PC123_meth_cpPCs_Neutro_genoPCs
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC1234_meth_cpPCs_Neutro_genoPCs.txt -Z -o KHS_sex_PC1234_meth_cpPCs_Neutro_genoPCs
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC12345_meth_cpPCs_Neutro_genoPCs.txt -Z -o KHS_sex_PC12345_meth_cpPCs_Neutro_genoPCs
/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k KHS_95_MAC02.BN.kinf -c covs_KHS_sex_PC13_meth_Neutro_genoPCs.txt -Z -o KHS_sex_PC13_meth_Neutro_genoPCs




/share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 -t beta_KHS -p age_KHS.txt -k  ninety_five_KHS_merged_52_tagged_real_segregating_MAF05_LD03.BN.kinf -c covs_KHS_sex_PC13_meth_Neutro_genoPCs.txt -Z -o KHS_sex_PC13_meth_Neutro_genoPCs









