
#  # in conda envt R_4.1.3

# #works on KHS subset, right now only BF_sig per CpG

for (pop in c("Himba")){
print(pop)
print(paste0("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/EMMAX_", pop, "_meths_updated_snp_ids_LD_pruned/95_MAC02_MeQTL_", pop, "_results"))
setwd(paste0("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/EMMAX_", pop, "_meths_updated_snp_ids_LD_pruned/95_MAC02_MeQTL_", pop, "_results"))

system("find . -type f -empty -exec rm {} +", intern = TRUE)

file_names <- list.files(pattern = "*.ps")

print(length(file_names))

# Create an empty data frame to store the results
result_df <- data.frame(cpg = character(), SNP = character(), effect = numeric (), P_value = numeric(), stringsAsFactors = FALSE)

# Iterate through the list of file names and read each file
for (file_name in file_names) {
if(file.exists(file_name)){
# Construct the full path to the file
full_path <- file.path(getwd(), file_name)

# Read the contents of the current file
file_content <- read.table(full_path, header = FALSE, stringsAsFactors = FALSE)
  
print(nrow((file_content)))

# Set the significance threshold
num_cis_snps <- nrow(file_content)
bf_sig <- .05/num_cis_snps 
bf_sig <- 1
#print(paste0("bf sig:", bf_sig))

# Check each row to see if the third column value is less than bf_sig
 significant_rows <- file_content[file_content$V3 < bf_sig, ]
 print(paste0("normal meqtl num sig: ", nrow(significant_rows)))
 # Add the significant rows to the result data frame with the file name in the first column
 if (nrow(significant_rows) > 0) {
   significant_rows$cpgs <- file_name
   #colnames(significant_rows)[1] <- "cpgs"
   print(significant_rows)
   result_df <- rbind(result_df, significant_rows[, c("cpgs", "V1", "V2","V3")])
 }
}
}
   
print(head(result_df))
write.table(result_df, paste0("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/", pop, "meqtls_all_tested.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep=",")}





