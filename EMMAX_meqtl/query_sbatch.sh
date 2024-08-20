#!/bin/bash
#SBATCH --job-name=all_process_chunks
#SBATCH --output=process_chunks_%A_%a.out
#SBATCH --error=process_chunks_%A_%a.err
#SBATCH --array=1-201000:1000  # Adjust based on the number of chunks
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=10:00:00
#SBATCH --mail-user=glmeeks@ucdavis.edu
#SBATCH --partition=production


# Population name
pop_name=$1

# Expected number of lines
expected_lines=$2

# Calculate the start and end line numbers for the current chunk
start=$((SLURM_ARRAY_TASK_ID))
end=$((start + 999))

start=1
end=$(wc -l < "${pop_name}_meqtls.txt")

echo "$start"
echo "$end"

# Process the current chunk
awk -v start="$start" -v end="$end" 'NR >= start && NR <= end' "${pop_name}_meqtls.txt" | awk -F ',' '{print $2}' | while IFS= read -r snp_name; do
    # Check if the file exists
    file_path="${pop_name}_meqtl_genos/${snp_name}_meqtl.txt"
    if [ ! -f "$file_path" ]; then
        echo "Processing SNP: $snp_name"
        echo "$file_path"
        grep "^$snp_name|" "${pop_name}_all_genos" > "$file_path"
    else
        # Check if the file has the expected number of lines
        num_lines=$(wc -l < "$file_path")
        if [ "$num_lines" -ne "$expected_lines" ]; then
            echo "Number of lines does not match the expected: $file_path"
            echo "Processing SNP: $snp_name"
            echo "$file_path"
            grep "^$snp_name|" "${pop_name}_all_genos" > "$file_path"
        else
            echo "File exists and has the expected number of lines: $file_path"
        fi
    fi
done


mv process_chunks_%A_%a.out slurm/
mv process_chunks_%A_%a.err slurm/
