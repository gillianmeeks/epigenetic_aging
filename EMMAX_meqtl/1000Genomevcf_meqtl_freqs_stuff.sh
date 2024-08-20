##make variant IDS without ref 
plink2 --vcf Baka_95_MAC02_pos_id.vcf --set-all-var-ids @:# -make-bed --out Baka_95_MAC02_pos_id_chrom_pos
plink2 --bfile Baka_95_MAC02_pos_id_chrom_pos --recode vcf --out Baka_95_MAC02_pos_id_chrom_pos 
plink2 --vcf Himba_95_MAC02_pos_id.vcf --set-all-var-ids @:# -make-bed --out Himba_95_MAC02_pos_id_chrom_pos
plink2 --bfile Himba_95_MAC02_pos_id_chrom_pos --recode vcf --out Himba_95_MAC02_pos_id_chrom_pos 
plink2 --vcf KHS_95_MAC02_rename_pos_ids.vcf --set-all-var-ids @:# -make-bed --out KHS_95_MAC02_rename_pos_ids_chrom_pos
plink2 --bfile KHS_95_MAC02_rename_pos_ids_chrom_pos --recode vcf --out KHS_95_MAC02_rename_pos_ids_chrom_pos 
plink2 --vcf 1kg_full_data_biallelic.vcf.gz  --set-all-var-ids @:# -make-bed --out 1kg_full_data_biallelic_chrom_pos
plink2 --bfile 1kg_full_data_biallelic_chrom_pos --recode vcf 1kg_full_data_biallelic_chrom_pos



#set same ref alt
plink2 --vcf KHS_95_MAC02_rename_pos_ids_chrom_pos.vcf --ref-from-fa /share/hennlab/users/glmeeks/age_methylation/methylation_imputation/saliva_imputation/imputation_processing/hs37d5.fa --make-bed -out KHS_95_MAC02_rename_pos_ids_chrom_pos_ref_from_fa
plink2 --bfile KHS_95_MAC02_rename_pos_ids_chrom_pos_ref_from_fa --recode vcf --out KHS_95_MAC02_rename_pos_ids_chrom_pos_ref_from_fa
plink2 --vcf Baka_95_MAC02_pos_id_chrom_pos.vcf --ref-from-fa /share/hennlab/users/glmeeks/age_methylation/methylation_imputation/saliva_imputation/imputation_processing/hs37d5.fa --make-bed -out Baka_95_MAC02_pos_id_chrom_pos_ref_from_fa
plink2 --bfile Baka_95_MAC02_pos_id_chrom_pos_ref_from_fa --recode vcf --out Baka_95_MAC02_pos_id_chrom_pos_ref_from_fa
plink2 --vcf Himba_95_MAC02_pos_id_chrom_pos.vcf --ref-from-fa /share/hennlab/users/glmeeks/age_methylation/methylation_imputation/saliva_imputation/imputation_processing/hs37d5.fa --make-bed -out Himba_95_MAC02_pos_id_chrom_pos_ref_from_fa
plink2 --bfile Himba_95_MAC02_pos_id_chrom_pos_ref_from_fa --recode vcf --out Himba_95_MAC02_pos_id_chrom_pos_ref_from_fa


##get positions for meqtls
cut -d "," -f 2 Baka_meqtls.txt | cut -d ":" -f 1 > Baka_meqtl_pos
cut -d "," -f 2 Baka_meqtls.txt | cut -d ":" -f 2 > Baka_meqtl_pos2
paste Baka_meqtl_pos Baka_meqtl_pos2 > Baka_meqtl_chr_pos
cut -d "," -f 2 KHS_meqtls.txt | cut -d ":" -f 1 > KHS_meqtl_pos
cut -d "," -f 2 KHS_meqtls.txt | cut -d ":" -f 2 > KHS_meqtl_pos2
paste KHS_meqtl_pos KHS_meqtl_pos2 > KHS_meqtl_chr_pos
cut -d "," -f 2 Himba_meqtls.txt | cut -d ":" -f 1 > Himba_meqtl_pos
cut -d "," -f 2 Himba_meqtls.txt | cut -d ":" -f 2 > Himba_meqtl_pos2
paste Himba_meqtl_pos Himba_meqtl_pos2 > Himba_meqtl_chr_pos

cat Baka_meqtl_chr_pos KHS_meqtl_chr_pos Himba_meqtl_chr_pos > pop_meqtls

#####subset to biallelic first
vcftools --gzvcf 1kg_full_data.vcf.gz --min-alleles 2 --max-alleles 2 --recode --out 1kg_full_data_biallelic
###restrict 1kg to those in the meqtl_chr_pos files ^
plink2 --vcf 1kg_full_data_biallelic.vcf --set-all-var-ids  @:# -make-bed -out 1kg_full_data_biallelic_chrom_pos
vcftools --vcf 1kg_full_data_biallelic.recode.vcf --positions pop_meqtls --recode --out 1kg_pop_meqtls

plink2 --vcf 1kg_pop_meqtls.recode.vcf  --set-all-var-ids  @:# -make-bed --out 1kg_pop_meqtls_chrom_pos
plink2 --bfile 1kg_pop_meqtls_chrom_pos --recode vcf --out 1kg_pop_meqtls_chrom_pos

plink2 --vcf 1kg_pop_meqtls_chrom_pos.vcf --ref-from-fa /share/hennlab/users/glmeeks/age_methylation/methylation_imputation/saliva_imputation/imputation_processing/hs37d5.fa --make-bed -out 1kg_pop_meqtls_chrom_pos_ref_from_fa
plink2 --bfile 1kg_pop_meqtls_chrom_pos_ref_from_fa --recode vcf --out 1kg_pop_meqtls_chrom_pos_ref_from_fa


##make the files with the unique snps that are meqtl i.e. so each snp is only listed once and extract these from the African pops
##extract freqs
vcftools --vcf  KHS_95_MAC02_rename_pos_ids_chrom_pos_ref_from_fa.vcf --positions unique_KHS_meqtl_chr_pos --recode --out KHS_95_MAC02_rename_pos_ids_chrom_pos_ref_from_fa_meqtls
plink2 --vcf KHS_95_MAC02_rename_pos_ids_chrom_pos_ref_from_fa_meqtls.recode.vcf --freq --out  KHS_95_MAC02_rename_pos_ids_chrom_pos_ref_from_fa_meqtls

plink --bfile Himba_95_MAC02_pos_id_chrom_pos_ref_from_fa --recode vcf --out Himba_95_MAC02_pos_id_chrom_pos_ref_from_fa
vcftools --vcf  Himba_95_MAC02_pos_id_chrom_pos_ref_from_fa.vcf --positions unique_Himba_meqtl_chr_pos --recode --out Himba_95_MAC02_pos_id_chrom_pos_ref_from_fa_meqtls
plink2 --vcf Himba_95_MAC02_pos_id_chrom_pos_ref_from_fa_meqtls.recode.vcf --freq --out  Himba_95_MAC02_pos_id_chrom_pos_ref_from_fa_meqtls

plink --bfile Baka_95_MAC02_pos_id_chrom_pos_ref_from_fa --recode vcf --out Baka_95_MAC02_pos_id_chrom_pos_ref_from_fa
vcftools --vcf  Baka_95_MAC02_pos_id_chrom_pos_ref_from_fa.vcf --positions unique_Baka_meqtl_chr_pos --recode --out Baka_95_MAC02_pos_id_chrom_pos_ref_from_fa_meqtls
plink2 --vcf Baka_95_MAC02_pos_id_chrom_pos_ref_from_fa_meqtls.recode.vcf --freq --out  Baka_95_MAC02_pos_id_chrom_pos_ref_from_fa_meqtls

plink2 --vcf 1kg_pop_meqtls_chrom_pos_ref_from_fa.recode.vcf --freq --out  1kg_pop_meqtls_chrom_pos_ref_from_fa_meqtls

####################FST###############
bgzip KHS_95_MAC02_rename_pos_ids_chrom_pos_ref_from_fa.vcf
bgzip Baka_95_MAC02_pos_id_chrom_pos_ref_from_fa.vcf
bgzip Himba_95_MAC02_pos_id_chrom_pos_ref_from_fa.vcf
bgzip 1kg_pop_meqtls_chrom_pos_ref_from_fa.vcf

bcftools index KHS_95_MAC02_rename_pos_ids_chrom_pos_ref_from_fa.vcf.gz
bcftools index Baka_95_MAC02_pos_id_chrom_pos_ref_from_fa.vcf.gz
bcftools index Himba_95_MAC02_pos_id_chrom_pos_ref_from_fa.vcf.gz
bcftools index 1kg_pop_meqtls_chrom_pos_ref_from_fa.vcf.gz

bcftools isec -p isec_KHS_1kg KHS_95_MAC02_rename_pos_ids_chrom_pos_ref_from_fa.vcf.gz 1kg_pop_meqtls_chrom_pos_ref_from_fa.vcf.gz
bcftools isec -p isec_Baka_1kg Baka_95_MAC02_pos_id_chrom_pos_ref_from_fa.vcf.gz 1kg_pop_meqtls_chrom_pos_ref_from_fa.vcf.gz
bcftools isec -p isec_Himba_1kg Himba_95_MAC02_pos_id_chrom_pos_ref_from_fa.vcf.gz 1kg_pop_meqtls_chrom_pos_ref_from_fa.vcf.gz

bgzip isec_KHS_1kg/0002.vcf 
bgzip isec_KHS_1kg/0003.vcf
bgzip isec_Baka_1kg/0002.vcf
bgzip isec_Baka_1kg/0003.vcf
bgzip isec_Himba_1kg/0002.vcf
bgzip isec_Himba_1kg/0003.vcf 
bcftools index isec_KHS_1kg/0002.vcf.gz
bcftools index isec_KHS_1kg/0003.vcf.gz
bcftools index isec_Baka_1kg/0002.vcf.gz
bcftools index isec_Baka_1kg/0003.vcf.gz
bcftools index isec_Himba_1kg/0002.vcf.gz
bcftools index isec_Himba_1kg/0003.vcf.gz

#create common variant vcfs between the pops and 1kg
bcftools merge isec_KHS_1kg/0002.vcf.gz isec_KHS_1kg/0003.vcf.gz -o common_KHS_1kg.vcf.gz
bcftools merge isec_Baka_1kg/0002.vcf.gz isec_Baka_1kg/0003.vcf.gz -o common_Baka_1kg.vcf.gz
bcftools merge isec_Himba_1kg/0002.vcf.gz isec_Himba_1kg/0003.vcf.gz -o common_Himba_1kg.vcf.gz

###fst and frequency example in himba
plink --vcf  common_Baka_1kg.vcf.gz --make-bed --out  common_Baka_1kg
paste Baka_inds Baka_inds | awk '{print $1, $2, "Baka"}' > Baka_plink_inds.txt
paste 1kg_inds 1kg_inds | awk '{print $1, $2, "1kg"}' > 1kg_plink_inds.txt
cat Baka_plink_inds.txt 1kg_plink_inds.txt > Baka_1kg_plink.txt
plink --vcf common_Baka_1kg --fst --within Baka_1kg_plink.txt --out Fst_Baka_1kg_meqtl













