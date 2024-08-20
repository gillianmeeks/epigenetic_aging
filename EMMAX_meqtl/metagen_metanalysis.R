
###in meta_R conda env in hennlab conda
##All 3 meta ##########  change to regress or no regress
library(meta)
#####################
args <- commandArgs(trailingOnly = TRUE)
j <- as.numeric(args[1])


load(paste0("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/fake_perms/Himba",j, "best_covs_with_SE_fake_meqtl_regressed.RData"))
#load("Himbabest_covs_with_SE_meqtl_regressed.RData")
#load("Himbabest_covs_with_SE_fusion_regressed.RData")
Himba_best_covs <- best_covs

load(paste0("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/fake_perms/Baka",j, "best_covs_with_SE_fake_meqtl_regressed.RData"))
#load("Bakabest_covs_with_SE_meqtl_regressed.RData")
#load("Bakabest_covs_with_SE_fusion_regressed.RData")
Baka_best_covs <- best_covs

load(paste0("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/fake_perms/KHS",j, "best_covs_with_SE_fake_meqtl_regressed.RData"))
#load("KHSbest_covs_with_SE_meqtl_regressed.RData")
#load("KHSbest_covs_with_SE_fusion_regressed.RData")
KHS_best_covs <- best_covs


Himba_best_covs$CPG.Labels <- as.character(Himba_best_covs$CPG.Labels)
KHS_best_covs$CPG.Labels <- as.character(KHS_best_covs$CPG.Labels)
Baka_best_covs$CPG.Labels <- as.character(Baka_best_covs$CPG.Labels)
common_sites <- intersect(intersect(Baka_best_covs$CPG.Labels, Himba_best_covs$CPG.Labels), KHS_best_covs$CPG.Labels)
#common_sites <- head(common_sites, 100)

#####restrict to initial meta sig sites#####
# load("metagen_sig_fixed_no_meqtl_regress.RData")
# meqtl_no_regress <- sig
# common_sites <- head(meqtl_no_regress$CPG.Labels, 6)
# print(common_sites)

###m.gen's pval.Q is the Test of heterogeniety p-value
#all 3
meta_gen <- data.frame()
for (i in common_sites){
    meta = data.frame()
    studlab <- c("KHS", "Baka", "Himba")
    TE<- c(KHS_best_covs[KHS_best_covs$CPG.Labels==i,"V2"], Baka_best_covs[Baka_best_covs$CPG.Labels==i, "V2"], Himba_best_covs[Himba_best_covs$CPG.Labels==i,"V2"])
    seTE <- c(KHS_best_covs[KHS_best_covs$CPG.Labels==i,"SE"], Baka_best_covs[Baka_best_covs$CPG.Labels==i, "SE"], Himba_best_covs[Himba_best_covs$CPG.Labels==i,"SE"])
    m.gen <- metagen(TE=TE, seTE=seTE, studlab=studlab, sm = "MD", fixed=TRUE, random=TRUE, method.tau ="REML", hakn=FALSE, prediction=TRUE, control=list(stepadj=0.5, maxiter=10000))
    new_row <- data.frame(CPG.Labels = i, P.value.random = m.gen$pval.random, P.value.fix= m.gen$pval.fixed, het.P.value = m.gen$pval.Q)
    meta_gen <- rbind(meta_gen, new_row)
    
}
print(nrow(meta_gen))
print(head(meta_gen))
save(meta_gen, file=paste0("fake_perms/metagen_all3_fake_meqtl_regressed", j, ".RData"))
#save(meta_gen, file="metagen_all3_meqtl_regressed.RData")
#save(meta_gen, file="metagen_all3_fusion_regressed.RData")

# #KHS and Baka
# meta_gen <- data.frame()
# for (i in common_sites){
#     meta = data.frame()
#     studlab <- c("KHS", "Baka")
#     TE<- c(KHS_best_covs[KHS_best_covs$CPG.Labels==i,"V2"], Baka_best_covs[Baka_best_covs$CPG.Labels==i, "V2"])
#     seTE <- c(KHS_best_covs[KHS_best_covs$CPG.Labels==i,"SE"], Baka_best_covs[Baka_best_covs$CPG.Labels==i, "SE"])
#     m.gen <- metagen(TE=TE, seTE=seTE, studlab=studlab, sm = "MD", fixed=TRUE, random=TRUE, method.tau ="REML", hakn=FALSE, prediction=TRUE, control=list(stepadj=0.5, maxiter=10000))
#     new_row <- data.frame(CPG.Labels = i, P.value.random = m.gen$pval.random, P.value.fix= m.gen$pval.fixed, het.P.value = m.gen$pval.Q)
#     meta_gen <- rbind(meta_gen, new_row)
    
# }
# print(nrow(meta_gen))
# save(meta_gen, file="metagen_baka_KHS.RData")

# #KHS and Himba
# meta_gen <- data.frame()
# for (i in common_sites){
#     meta = data.frame()
#     studlab <- c("KHS", "Himba")
#     TE<- c(KHS_best_covs[KHS_best_covs$CPG.Labels==i,"V2"], Himba_best_covs[Himba_best_covs$CPG.Labels==i,"V2"])
#     seTE <- c(KHS_best_covs[KHS_best_covs$CPG.Labels==i,"SE"], Himba_best_covs[Himba_best_covs$CPG.Labels==i,"SE"])
#     m.gen <- metagen(TE=TE, seTE=seTE, studlab=studlab, sm = "MD", fixed=TRUE, random=TRUE, method.tau ="REML", hakn=FALSE, prediction=TRUE, control=list(stepadj=0.5, maxiter=10000))
#     new_row <- data.frame(CPG.Labels = i, P.value.random = m.gen$pval.random, P.value.fix= m.gen$pval.fixed, het.P.value = m.gen$pval.Q)
#     meta_gen <- rbind(meta_gen, new_row)
    
# }
# print(nrow(meta_gen))
# save(meta_gen, file="metagen_KHS_himba.RData")

# #baka and himba
# meta_gen <- data.frame()
# for (i in common_sites){
#     meta = data.frame()
#     studlab <- c("Baka", "Himba")
#     TE<- c(Baka_best_covs[Baka_best_covs$CPG.Labels==i, "V2"], Himba_best_covs[Himba_best_covs$CPG.Labels==i,"V2"])
#     seTE <- c(Baka_best_covs[Baka_best_covs$CPG.Labels==i, "SE"], Himba_best_covs[Himba_best_covs$CPG.Labels==i,"SE"])
#     m.gen <- metagen(TE=TE, seTE=seTE, studlab=studlab, sm = "MD", fixed=TRUE, random=TRUE, method.tau ="REML", hakn=FALSE, prediction=TRUE, control=list(stepadj=0.5, maxiter=10000))
#     new_row <- data.frame(CPG.Labels = i, P.value.random = m.gen$pval.random, P.value.fix= m.gen$pval.fixed, het.P.value = m.gen$pval.Q)
#     meta_gen <- rbind(meta_gen, new_row)
    
# }
# print(nrow(meta_gen))
# save(meta_gen, file="metagen_a.RData")






















