load("/share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/Diverse_Age_QC_methylation_phenos.RData")
#load("/share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/EMMAX/metagen_sig_fixed.RData") 
library("glmnet")
manifest_450k <- read.csv("/share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/humanmethylation450_15017482_v1-2.csv", sep=",", header=TRUE)
manifest_EPIC <- read.csv("/share/hennlab/data/methylation_data/Himba/EPIC32_042020/ReferenceFiles/MethylationEPIC_v-1-0_B4_noheading.csv", header=TRUE,sep=",")

#library(ggplot2)
library(data.table)
library(dplyr)

process_column_name <- function(cpg) {
  cpg_ <- strsplit(cpg, split = "_")[[1]][2]
  cpg_ <- strsplit(cpg_, split = "\\.txt")[[1]][1]  
  return(cpg_)
}


args <- commandArgs(trailingOnly = TRUE)
j <- as.numeric(args[1])
pop <- as.character(args[2])
print(pop)

# Load phenotype data based on the population
if (pop == "KHS") {
  pheno <- KHS_pheno_merged
  rownames(pheno) <- pheno$IID
  manifest <- manifest_450k 
  tfam <- read.delim("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/beta_meqtl_regress_KHS.tfam", sep='', header=FALSE)
  n <- 52
} else if (pop == "Baka") {
  pheno <- Baka_pheno
  rownames(pheno) <- pheno$IID
  manifest <- manifest_450k 
    tfam <- read.delim("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/beta_meqtl_regress_Baka.tfam", sep='', header=FALSE)
  n <- 35
} else if (pop == "Himba") {
  pheno <- Himba_pheno_merged
  rownames(pheno) <- pheno$ID
  pheno$IID <- pheno$ID
  pheno <- pheno[!pheno$ID == "HMB181_2", ]
  manifest <- manifest_EPIC
  n <- 51
 tfam <- read.delim("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/beta_meqtl_regress_Himba.tfam" ,sep='', header=FALSE)
    
} else if (pop == "all") {
  load("/share/hennlab/users/glmeeks/age_methylation/age_prediction/merged_data_all_phenos.RData")
  pheno <- pheno_merged
  rownames(pheno) <- pheno$ID
  pheno$IID <- pheno$ID
  manifest <- manifest_450k
  n <- 138
  tfam <-read.delim("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/beta_all_no_regress.tfam", sep='', header=FALSE)
}

# Set working directory
setwd("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/")

# Read methylation data
meth_data <- fread(paste0(pop, "_meth_phenos.tsv"))  
row_index <- which(colnames(meth_data) == "HMB494.2")
names(meth_data)[row_index] <- "HMB494-2"
pheno <- pheno[c(colnames(meth_data[,5:ncol(meth_data)])),]
rownames(meth_data) <- meth_data$cpg
residuals_df <- t(meth_data[,5:ncol(meth_data)]) 
colnames(residuals_df) <- meth_data$cpg

meqtls_df <- read.table(paste0(pop, "_meqtls.txt"), sep = ",", header = FALSE)
meqtls_df$V1 <- sapply(meqtls_df[,1], process_column_name)
meqtls_df <- meqtls_df[!duplicated(meqtls_df$V1), ]

meqtls_df$fake_meqtl <- sapply(1:nrow(meqtls_df), function(i) {
  curr_chrom <- meqtls_df[i, "chrom"]
  available_values <- meqtls_df[meqtls_df$chrom != curr_chrom, "V2"]
  sampled_value <- if(length(available_values) > 0) sample(available_values, 1) else NA
  return(sampled_value)
})

meqtls_df$chrom <- sapply(meqtls_df$V2, function(x) {
  unlist(strsplit(x, split = ":"))[1]
})

unique_chroms <- unique(meqtls_df$chrom)
sampled_values <- numeric(nrow(meqtls_df))
for (chrom in unique_chroms) {
  idx <- meqtls_df$chrom != chrom
  sampled_values[idx] <- sample(meqtls_df$V2[idx], sum(idx), replace = TRUE)
}

meqtls_df$fake_meqtl <- sampled_values


perform_regression <- function(cpg_) {
    meth <- as.data.frame(residuals_df[, colnames(residuals_df) == cpg_])
    if (cpg_ %in% meqtls_df$V1) {
        print(cpg_)
        meth$ID <- rownames(meth)  
        #print(meqtls_df[meqtls_df$V1 == cpg_,])
        snp_name <- meqtls_df[meqtls_df$V1 == cpg_, "fake_meqtl"]
        file_path <- paste0(pop, "_meqtl_genos/", snp_name, "_meqtl.txt")
        snp_data <- read.delim(file_path, header = FALSE, sep = "|")
        snp_data$sum_alleles <- sapply(snp_data$V7, function(x) {
            if (x == "0/0") {
                return(0)
            } else if (x == "0/1") {
                return(1)
            } else if (x == "1/1") {
                return(2)
            } else {
                return(NA) 
            }
        })
        snp_data$V6 <- gsub(".*_", "", snp_data$V6)
        rownames(snp_data) <- snp_data$V6
        snp_data <- snp_data[rownames(meth), ]
        pheno <- pheno[rownames(meth), ]
        snp_data$ID <- rownames(snp_data)
        snp_data$meth_value <- meth[match(snp_data$ID, meth$ID), 1]
        snp_data$age <- pheno[match(snp_data$ID, pheno$IID), "age"]
        model_regress <- lm(meth_value ~ sum_alleles, data = snp_data, na.action = na.exclude)
        to_return <- as.data.frame(resid(model_regress))
        names(to_return) <- cpg_
        return(to_return)
    } else {
        to_return <- meth
        names(to_return) <- cpg_
        return(to_return)
    }
}

#############
##subsample##
#############
#residuals_df <- residuals_df[,1:100]
#############
# Apply perform_regression function to every column of residuals_df
residuals_df <- as.data.frame(lapply(colnames(residuals_df), function(cpg_) {
    perform_regression(cpg_)
}))


##convert to tped####
            tped  <- residuals_df  ##na.omit the first run so I should do EMMAX again with the no residuals
            ones <- rep(1, nrow(pheno))

            ####ORDER of inds in pheno and meth matrix must be same as geno tfam
          
            tped <- as.data.frame(t(tped))
            rownames(tped) <- colnames(residuals_df)
            print(head(tped))
            print(head(tfam))
            manifest <- manifest[, c("IlmnID", "CHR", "MAPINFO")]
            tped <- tped[rownames(tped)%in% manifest$IlmnID,]
            tped <- tped[,c(tfam$V1)]
            pheno <- pheno[c(tfam$V1),]  
            print(head(pheno[, 1:3]))
            ######
            siteinfo <- manifest[manifest$IlmnID %in% rownames(tped),] 
            tped$CpG.Labels <- rownames(tped)
            siteinfo <- dplyr::arrange(siteinfo, desc(IlmnID))
            tped<- dplyr::arrange(tped, desc(CpG.Labels))
            CpG.Labels <- rownames(tped)
            chr <- siteinfo[match(tped$CpG.Labels, siteinfo$IlmnID), "CHR"]
            gen_pos <- rep(0, nrow(tped))
            phy_pos <- siteinfo[match(tped$CpG.Labels, siteinfo$IlmnID), "MAPINFO"]
            tped <- tped[,!colnames(tped)%in% c("CpG.Labels")]
            tped <- tped[,c(tfam$V1)]
            tped <- cbind(chr, CpG.Labels, gen_pos, phy_pos, tped)
	    print(head(tped))       
write.table(tped, file = paste0("fake_perms/", j,"beta_fake_meqtl_regress_", pop, ".tped"), quote = FALSE, col.names = FALSE, row.names = FALSE)
            print("done tped")  
write.table(tfam, file=paste0("fake_perms/", j,"beta_fake_meqtl_regress_", pop, ".tfam"), quote=FALSE, col.names=FALSE, row.names=FALSE)
	    print(head(tfam)) 
	    print("done tfam")
