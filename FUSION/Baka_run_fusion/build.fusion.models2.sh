#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH --mem=30Gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=glmeeks@ucdavis.edu
#SBATCH  --array=1-400900:100

##SBATCH --array=1-1000:100

#SBATCH --partition=production
#SBATCH --job-name=Baka_fusion
# Define output and error filenames
output_file="${SLURM_ARRAY_TASK_ID}A.out"
error_file="${SLURM_ARRAY_TASK_ID}A.err"
#SBATCH --output="$output_file"
#SBATCH --error="$error_file"

###set 8 hours and 20Gb for himba, set 6 hours and 16Gb for others
##CANT RUN MULTIPLE AT ONCE IN THE SAME DIRECTORY
##not inclusive of the y in x-y
##sbatch build.fusion.models2.sh Himba #cpgs
##change job-name
source /share/hennlab/progs/miniconda3/etc/profile.d/conda.sh
conda activate FUSION_env
#Allow for NA cpgs
#713771, 710576 no NA cpgs for himba (Himba have 217 fewer probes when you overlay on manifest, I think because ctrl probes were still in)
#400893, 399543 cpgs for Baka
#418629, 416692 for KHS
#355103, 352053 in all 


pop="$1"
cpgs="$2"

# Check if Conda environment was activated successfully
if [ $? -eq 0 ]; then
    echo "Conda environment 'FUSION_env' activated successfully."
else
    echo "Error: Failed to activate Conda environment 'FUSION_env'."
    exit 1
fi



## BUILD WEIGHTS FOR ALL GENES WITH SIG. H2 IN EITHER EA OR AA

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

GCTA=/software/gcta/1.93.1beta/lssc0-linux/gcta
GEMMA=/software/gemma/0.98.1/lssc0-linux/bin/gemma 
PLINK=/software/plink/1.90p/x86_64-linux-ubuntu14.04/bin/plink

# Specify directory and data locations
DATA_DIR=/share/hennlab/users/glmeeks/age_methylation/FUSION/FUSION_data/
COV_SUBDIR=covar/
EXPR_SUBDIR=FUSION_pheno/
OUT_DIR=/share/hennlab/users/glmeeks/age_methylation/FUSION/FUSION_final_SAVE_HSQ/
RUN_DIR=/share/hennlab/users/glmeeks/age_methylation/FUSION/${pop}_run_fusion/

FUSION=${RUN_DIR}/compute_weights_wSusie.R

# Specify usage limits for PLINK
MEM=4000
N_THREADS=1


start=$NR
stop=$((NR + 99))
P=$1
 
cd $RUN_DIR
for IDX in `seq $start $stop`
do
    # Extract parameters for FUSION run
    
    params=`sed "${IDX}q;d" ${RUN_DIR}${P}_fusion_params2.tsv`
    #params=`sed "${IDX}q;d" ${RUN_DIR}not_done_${P}_fusion_params2.tsv`
        set -- junk $params
        shift
    CHR=$2
    P0=`echo $params | awk '{ p=$3 - 500e3; if(p<0) p=0; print p; }'`
    P1=`python -c "print( int(int($4) + 500e3))"`
    GENE=$5
    echo SUMMARY ${P}, ${GENE}, ${CHR}, ${P0}, ${P1}

    # Point to genotype data
    BFILE=${DATA_DIR}/geno/${P}/${P}_chr${CHR}_FUSION_data
    # Point to covariate files for this tissue
    COVARS=${DATA_DIR}${COV_SUBDIR}/${P}_covars.tsv
     
    echo ${COVARS}
    echo ${BFILE}
    # tempdir for model generation
    TMPDIR=temp_$IDX
    # remove if exists, then create
    if [ -d $TMPDIR ]; then
        rm -fr $TMPDIR
    fi
    mkdir -p $TMPDIR
    # Define temporary outfile names
    OUT=$TMPDIR/${IDX}

    # START TIMING
    start_time=$SECONDS

    # Create a temporary header file that contains the sample IDs
    zcat ${DATA_DIR}/${EXPR_SUBDIR}/${P}/${P}_hg19_meth.bed.gz | head -n1 | tr '\t' '\n' > $TMPDIR/header
    #zcat ${DATA_DIR}/${EXPR_SUBDIR}/${P}/with_header_not_done_${P}_hg19_meth.bed.gz | head -n1 | tr '\t' '\n' > $TMPDIR/header


    # Create phenotype file for PLINK
    #expr=`zgrep -w $GENE ${DATA_DIR}/${EXPR_SUBDIR}/${P}/with_header_not_done_${P}_hg19_meth.bed.gz`
    expr=`zgrep -w $GENE ${DATA_DIR}/${EXPR_SUBDIR}/${P}/${P}_hg19_meth.bed.gz`
    echo $expr | tr ' ' '\n' | paste $TMPDIR/header $TMPDIR/header - | tail -n+5 > ${OUT}.pheno
    echo `head ${OUT}.pheno`

    # Extract genomic region
    $PLINK --allow-no-sex --keep-allele-order --bfile $BFILE --chr $CHR --from-bp $P0 --to-bp $P1 --maf 0.0001 \
        --snps-only --keep $OUT.pheno --pheno $OUT.pheno --threads $N_THREADS --memory $MEM --make-bed --out $OUT \
        --silent

    # if doesn't exist something went wrong; continue
    if [ ! -f $OUT.bed ]
    then
        continue
        echo "no plink output"
    fi
    cd $TMPDIR
    OUT=$IDX
    ln -s ./ output

    # Make directory to store results
    echo ""${OUT_DIR}/${P}/${GENE}""
    mkdir -p ${OUT_DIR}/${P}/${GENE}

    # Run FUSION 
    Rscript $FUSION \
        --bfile $OUT \
        --pheno $OUT.pheno \
        --tmp TEMP \
        --out ${OUT_DIR}/${P}/${GENE}/${GENE} \
        --save_hsq TRUE \
        --hsq_p 0.05 \
        --PATH_gcta $GCTA \
        --PATH_gemma $GEMMA \
        --PATH_plink $PLINK \
        --plink_threads $N_THREADS \
        --plink_mem $MEM \
        --models enet,susie,top1,lasso \
        --covar $COVARS \
        --verbose 2
   echo  "done with FUSION $FUSION $OUT --out ${OUT_DIR}/${P}/${GENE}/${GENE} --covar $COVARS" 
   

    # Move back up to run directory, clean up temporary directory
    cd $RUN_DIR
    rm -f -r $TMPDIR
    end=$SECONDS
    TIME_SEC=`echo $end - $start_time | bc | cat`
    echo "SUMMARY |  $P, $GENE, $CHR, $P0, $P1, $TIME_SEC s total time"
done


#job_name="${SLURM_JOB_NAME}"

#rm "resources_"${SLURM_ARRAY_TASK_ID}".txt"
#touch "resources_"${SLURM_ARRAY_TASK_ID}".txt"



#sacct -j ${SLURM_ARRAY_TASK_ID} --format=JobID,JobName,Elapsed,MAXRSS,ExitCode>> "resources_"${SLURM_ARRAY_TASK_ID}".txt"


