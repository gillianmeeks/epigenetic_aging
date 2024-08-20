#!/bin/sh
##TEST git####
##allowed for nas here
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
##SBATCH --time=8:00:00
#SBATCH --mem=10Gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=glmeeks@ucdavis.edu

##SBATCH --array=1-401000:1000
##SBATCH --array=1-3000:1000
##SBATCH --job-name=Baka_EMMAX

##SBATCH --array=1-714000:1000
##SBATCH --array=1-3000:1000
##SBATCH --job-name=Himba_EMMAX

#SBATCH --array=1-418700:1000
##SBATCH --array=1-3000:1000
#SBATCH --job-name=KHS_EMMAX

##SBATCH  --array=1-355200:1000
##SBATCH --job-name=all_EMMAX_rest

#SBATCH --partition=production


PLINK=/software/plink/1.90p/x86_64-linux-ubuntu14.04/bin/plink

##ran this including NA cpgs
##not inclusive of the y in x-y
##change job-name
#713771 (710576 no NAs) cpgs for himba (Himba have 217 fewer probes when you overlay on manifest, not sure why these are missing from the manifest)
#400893 (399543 no NAs) cpgs for Baka
#418629 (416692 no NAs) for KHS
#355103,(352053 no NAs) in all 

mem_before=$(ps -p $$ -o rss=)
start_time=$(date +%s)


pop="$1"

outfile="${SLURM_ARRAY_TASK_ID}A.out"
error_file="${SLURM_ARRAY_TASK_ID}A.err"
#SBATCH --output="$outfile"
#SBATCH --error="$error_file"


cpgs=`wc -l < ""$pop"_meth_phenos.tsv"`
#cpgs=`wc -l < "with_header_not_done_${pop}_meth_phenos.tsv"`

echo ""$cpgs""

PLINK=/software/plink/1.90p/x86_64-linux-ubuntu14.04/bin/plink
#LD pruned file at 50 snp, moving 5, r2 .5
#BFILE=/share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/EMMAX/${pop}_95_MAC02_ids_LD_50_5_50   

#fixed rsids to all be pos:ref:alt
BFILE=/share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/EMMAX/${pop}_95_MAC02_pos_ids_LD_50_5_5

#BFILE=/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/all_pops_merged_segregating_pos_ids_maf_01_LD_50_5_05

if [ ! $SLURM_ARRAY_TASK_ID ]; then
    NR=$1
    SCRATCHDIR=`echo $PWD`
else
    NR=$SLURM_ARRAY_TASK_ID
fi
#this is where you can set what line to go to
if [ $NR -gt $cpgs ]; then
    break
fi

echo "task ID: $SLURM_ARRAY_TASK_ID"
echo "line: $NR"

start=$NR
stop=$((NR + 999))

mkdir -p "EMMAX_"$pop"_meths_updated_snp_ids_LD_pruned"
cd "EMMAX_"$pop"_meths_updated_snp_ids_LD_pruned"
#mkdir -p "Full_95_MAC02_MeQTL_"$pop"_results" 
mkdir -p "95_MAC02_MeQTL_"$pop"_results" 

# Extract the first row from the ${pop}_meth_phenos.tsv and store it in a variable

header=$(head -1 "../${pop}_meth_phenos.tsv")


# Extract the first two columns from the header (excluding the first four column names)
echo "$header" | awk '{$1=$2=$3=$4=""; print $0}' | tr -s ' '  '\n' > "temp${SLURM_ARRAY_TASK_ID}"
paste "temp${SLURM_ARRAY_TASK_ID}" "temp${SLURM_ARRAY_TASK_ID}" > "temp_${SLURM_ARRAY_TASK_ID}"
cat "temp_${SLURM_ARRAY_TASK_ID}" | head 

for line_number in $(seq $start $stop); do
    echo "${line_number} line number"   
    # Get the corresponding cpg line from the ${pop}_meth_phenos.tsv
    cpg=$(tail -n +2 "../${pop}_meth_phenos.tsv" | sed -n "${line_number}p" )
    #cpg=$(tail -n +2 "../${pop}_meth_phenos.tsv" | sed -n "${line_number}p" )
    CHR=$(echo "$cpg" | cut -f 1)
    P0=$(echo "$cpg" | awk '{ p=$2 - 100e3; if(p<0) p=0; print p; }')
    P1=$(echo "$cpg" | awk '{ p=$3 + 100e3; if(p<0) p=0; print p; }')
    CPG=$(echo "$cpg" | cut -f 4)
    echo SUMMARY ${CPG}, ${CHR}, ${P0}, ${P1}
    output_file="cpg_${CPG}.txt"
  
    # Append the cpg line (excluding the first four columns) to the output file
    echo "$cpg" | awk '{$1=$2=$3=$4=""; print $0}' | tr -s ' '  '\n' > "pheno${line_number}"
    paste "temp_${SLURM_ARRAY_TASK_ID}" "pheno${line_number}" > "space${line_number}" 
    rm "pheno${line_number}"
    echo "number phenos $(wc -l <  "space${line_number}")"
    tail -n +2  "space${line_number}" > "$output_file"
    cat "cpg_${CPG}.txt" | head
    rm "space${line_number}" 
   # Extract genomic region
    cat $output_file | cut -f 1-2 > "KEEP${line_number}"
    cat "KEEP${line_number}" | head
    echo "number phenos $(wc -l < "KEEP${line_number}")"

    echo "    
     "$PLINK" --allow-no-sex --keep-allele-order --bfile "$BFILE" --chr "$CHR" --from-bp "$P0" --to-bp "$P1" \
         --snps-only --keep KEEP"$line_number" --silent --make-bed --out ""$pop"_"$CPG"_cis""
    
     $PLINK --allow-no-sex --keep-allele-order --bfile $BFILE --chr $CHR --from-bp $P0 --to-bp $P1 \
         --snps-only --keep "KEEP${line_number}" --silent --make-bed --out ""$pop"_"$CPG"_cis" 
        
    rm "KEEP${line_number}"
    "$PLINK" --bfile ""$pop"_"$CPG"_cis" --recode 12 --output-missing-genotype 0 --transpose --out ""$pop"_"$CPG"_cis"
    
    
    /share/hennlab/progs/emmax-beta-07Mar2010/emmax -v -d 10 \
    -t ""$pop"_"$CPG"_cis" \
    -p "$output_file" \
    -k "/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/"$pop"_95_MAC02_pos_id.BN.kinf" \
    -c "../covs_"$pop"_with_age.txt" \
    -o "95_MAC02_MeQTL_"$pop"_results/"$output_file"_meqtl" 
    
    rm  ""$pop"_"$CPG"_cis.bed"
    rm  ""$pop"_"$CPG"_cis.bim"
    rm  ""$pop"_"$CPG"_cis.fam"
    rm  ""$pop"_"$CPG"_cis.tfam"
    rm  ""$pop"_"$CPG"_cis.tped"
    rm  ""$pop"_"$CPG"_cis.log"
    rm ""$pop"_"$CPG"_cis.nosex"
    rm "cpg_${CPG}.txt"
    rm   "95_MAC02_MeQTL_"$pop"_results/"$output_file"_meqtl.reml"
    rm   "95_MAC02_MeQTL_"$pop"_results/"$output_file"_meqtl.log"


done
rm "temp${SLURM_ARRAY_TASK_ID}"
rm "temp_${SLURM_ARRAY_TASK_ID}"
# Calculate the end time after the loop finishes
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$((elapsed_time % 3600 / 60))
seconds=$((elapsed_time % 60))

# Calculate the memory usage after the operations
mem_after=$(ps -p $$ -o rss=)
mem_diff=$((mem_after - mem_before))

# Echo the total elapsed time and memory usage to the .out file
echo "Total Elapsed Time: ${hours}h:${minutes}m:${seconds}s" 
echo "Total Memory Used: ${mem_diff} KB"



#after this run make_meqtl.R with changed pop argument to get the csv of meqtls  and then query.sh to get the genoyptes for these meqtls.
