load("/share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/Diverse_Age_QC_methylation_phenos.RData")
#load("/share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/EMMAX/metagen_sig_fixed.RData") 
library("glmnet")
library("dplyr")
manifest_450k <- read.csv("/share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/humanmethylation450_15017482_v1-2.csv", sep=",", header=TRUE)
manifest_EPIC <- read.csv("/share/hennlab/data/methylation_data/Himba/EPIC32_042020/ReferenceFiles/MethylationEPIC_v-1-0_B4_noheading.csv", header=TRUE,sep=",")

library(data.table)
library(dplyr)

process_column_name <- function(cpg) {
  cpg_ <- strsplit(cpg, split = "_")[[1]][2]
  cpg_ <- strsplit(cpg_, split = "\\.txt")[[1]][1]  
  return(cpg_)
}


args <- commandArgs(trailingOnly = TRUE)
pop <- as.character(args[1])
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
no_residuals_df <- t(meth_data[,5:ncol(meth_data)]) 
colnames(residuals_df) <- meth_data$cpg
colnames(no_residuals_df) <- meth_data$cpg

meqtls_df <- read.table(paste0(pop, "_meqtls.txt"), sep = ",", header = FALSE)
meqtls_df$V1 <- sapply(meqtls_df[,1], process_column_name)

###only regress out lowest p-value snp##
meqtls_df <- meqtls_df %>%
  group_by(V1) %>%
  filter(V4 == min(V4)) %>%
  slice(1) %>%
  ungroup()
 
write.table(meqtls_df, file=paste0("top_meqtl_", pop, ".txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)

# perform_regression <- function(cpg_) {
#     meth <- as.data.frame(residuals_df[, colnames(residuals_df) == cpg_])
#     if (cpg_ %in% meqtls_df$V1) {
#         print(cpg_)
#         meth$ID <- rownames(meth)  
#         print(meqtls_df[meqtls_df$V1 == cpg_,])
#         snp_name <- meqtls_df[meqtls_df$V1 == cpg_, 2]
#         file_path <- paste0(pop, "_meqtl_genos/", snp_name, "_meqtl.txt")
#         snp_data <- read.delim(file_path, header = FALSE, sep = "|")
#         snp_data$sum_alleles <- sapply(snp_data$V7, function(x) {
#             if (x == "0/0") {
#                 return(0)
#             } else if (x == "0/1") {
#                 return(1) 
#             } else if (x == "1/1") {
#                 return(2)
#             } else {
#                 return(NA) 
#             }
#         })
#         snp_data$V6 <- gsub(".*_", "", snp_data$V6)
#         rownames(snp_data) <- snp_data$V6
#         snp_data <- snp_data[rownames(meth), ]
#         pheno <- pheno[rownames(meth), ]
#         snp_data$ID <- rownames(snp_data)
#         snp_data$meth_value <- meth[match(snp_data$ID, meth$ID), 1]
#         snp_data$age <- pheno[match(snp_data$ID, pheno$IID), "age"]
#         model_regress <- lm(meth_value ~ sum_alleles, data = snp_data, na.action = na.exclude)
#         to_return <- as.data.frame(resid(model_regress))
#         names(to_return) <- cpg_
#         return(to_return)
#     } else {
#         print(cpg_)
#         print("no meqtl")
#         to_return <- meth
#         names(to_return) <- cpg_
#         return(to_return)
#     }
# }  

# #############
# ##subsample##
# #############
# #residuals_df <- residuals_df[,1:1000]
# #############

# residuals_df <- as.data.frame(lapply(colnames(residuals_df), function(cpg_) {
#     perform_regression(cpg_)
# }))
#     save(residuals_df, file=paste0("/share/hennlab/users/glmeeks/age_methylation/age_prediction/", pop, "_top_meqtl_all_sites_regressed.RData"))
#     save(no_residuals_df, file=paste0("/share/hennlab/users/glmeeks/age_methylation/age_prediction/", pop, "_all_sites_no_regressed.RData"))
   
#     ###convert to tped####
#     tped  <- residuals_df  
#     tped2 <- no_residuals_df
#     ones <- rep(1, nrow(pheno))
#     tped <- as.data.frame(t(tped))
#     tped2 <- as.data.frame(t(tped2))
#     rownames(tped) <- colnames(residuals_df)
#     rownames(tped2) <- colnames(no_residuals_df)
#     tped <- tped[,c(tfam$V1)]
#     tped2 <- tped2[,c(tfam$V1)]
#     pheno <- pheno[c(tfam$V1),]  
#     manifest <- manifest[, c("IlmnID", "CHR", "MAPINFO")]
#     tped <- tped[rownames(tped) %in% manifest$IlmnID,]
#     tped2 <- tped2[rownames(tped2) %in% manifest$IlmnID,]
#     print(nrow(tped))
#     siteinfo <- manifest[manifest$IlmnID %in% rownames(tped),] 
#     siteinfo2 <- manifest[manifest$IlmnID %in% rownames(tped2),]   
#     tped$CpG.Labels <- rownames(tped)
#     tped2$CpG.Labels <- rownames(tped2)
#     siteinfo <- dplyr::arrange(siteinfo, desc(IlmnID))
#     siteinfo2 <- dplyr::arrange(siteinfo2, desc(IlmnID))
#     tped<- dplyr::arrange(tped, desc(CpG.Labels))
#     tped2<- dplyr::arrange(tped2, desc(CpG.Labels))
#     print(identical(tped$CpG.Labels, siteinfo$IlmnID))
#     print(identical(tped2$CpG.Labels, siteinfo2$IlmnID))
#     CpG.Labels <- rownames(tped)
#     CpG.Labels <- rownames(tped2)
#     #head(siteinfo[,c("IlmnID", "CHR", "MAPINFO")])
#     chr <- siteinfo[match(tped$CpG.Labels, siteinfo$IlmnID), "CHR"]
#     gen_pos <- rep(0, nrow(tped))
#     phy_pos <- siteinfo[match(tped$CpG.Labels, siteinfo$IlmnID), "MAPINFO"]
#     tped <- tped[,!colnames(tped)%in% c("CpG.Labels"),]
#     tped <- tped[,c(tfam$V1)]
#     tped <- cbind(chr, CpG.Labels, gen_pos, phy_pos, tped)   
#     chr <- siteinfo2[match(tped2$CpG.Labels, siteinfo2$IlmnID), "CHR"]
#     gen_pos <- rep(0, nrow(tped2))
#     phy_pos <- siteinfo2[match(tped2$CpG.Labels, siteinfo2$IlmnID), "MAPINFO"]
#     tped2 <- tped2[,!colnames(tped2)%in% c("CpG.Labels")]
#     tped2 <- tped2[,c(tfam$V1)]
#     tped2 <- cbind(chr, CpG.Labels, gen_pos, phy_pos, tped2)
       
#     write.table(tped, file=paste0("testbeta_meqtl_regress_", pop, ".tped"), quote=FALSE, col.names=FALSE, row.names=FALSE)
#     write.table(tped2, file=paste0("testbeta_no_regress_", pop, ".tped"), quote=FALSE, col.names=FALSE, row.names=FALSE)
#     write.table(tfam, file=paste0("testbeta_meqtl_regress_", pop, ".tfam"), quote=FALSE, col.names=FALSE, row.names=FALSE)
#     write.table(tfam, file=paste0("testbeta_no_regress_", pop, ".tfam"), quote=FALSE, col.names=FALSE, row.names=FALSE)


