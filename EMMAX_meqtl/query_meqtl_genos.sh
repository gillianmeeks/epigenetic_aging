
#run this first to get all genos, fast step

#bcftools view /share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/EMMAX/Himba_95_MAC02_pos_ids.vcf.gz  | bcftools query -f'[%ID|%CHROM|%POS|%REF|%ALT|%SAMPLE|%GT\n]' > "Himba_all_genos"


#bcftools view /share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/EMMAX/KHS_95_MAC02_rename_pos_ids.vcf.gz | bcftools query -f'[%ID|%CHROM|%POS|%REF|%ALT|%SAMPLE|%GT\n]' > "KHS_all_genos"


#bcftools view /share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/EMMAX/Baka_95_MAC02_pos_ids.vcf.gz | bcftools query -f'[%ID|%CHROM|%POS|%REF|%ALT|%SAMPLE|%GT\n]' > "Baka_all_genos"

#bcftools view /share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/all_pops_merged_segregating_pos_ids_maf_01.vcf.gz | bcftools query -f'[%ID|%CHROM|%POS|%REF|%ALT|%SAMPLE|%GT\n]' > "all_all_genos"
pop="$1"

# then get only significant ones in R into a csv
# run this
head "${pop}_meqtls.txt"
file="${pop}_meqtls.txt"
mkdir -p "${pop}_meqtl_genos"

awk -F ',' '{print $2}' "$file" | while IFS= read -r snp_name; do
    # Check if the file exists
    if [ ! -e "${pop}_meqtl_genos/${snp_name}_meqtl.txt" ]; then
        echo "Processing SNP: $snp_name"
        grep "^"$snp_name"|" "${pop}_all_genos" > "${pop}_meqtl_genos/${snp_name}_meqtl.txt"
    else
        echo "File already exists: ${pop}_meqtl_genos/${snp_name}_meqtl.txt"
    fi
done

# pop="$1"

# #then get only sig ones in R into a csv, then run this 
# head ""$pop"_age_meqtls.txt"
# file=""$pop"_age_meqtls.txt"
# mkdir -p ""$pop"_meqtl_genos"


# awk -F ',' '{print $2}' $file | while IFS= read -r snp_name; do
#     echo "Processing SNP: $snp_name"
#     grep "^$snp_name" ""$pop"_all_genos" > ""$pop"_meqtl_genos/${snp_name}_meqtl.txt"
# done



