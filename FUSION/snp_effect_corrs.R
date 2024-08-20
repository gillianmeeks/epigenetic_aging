args <- commandArgs(trailingOnly = TRUE)
pop1 <- as.character(args[1])
pop2 <- as.character(args[2])
#pop1 <- "KHS"
#pop2 <- "Baka"

    population_data1 <- read.delim(paste0("FUSION_final_SAVE_HSQ/",pop1,"_R2_hsq_4.tsv"))
    population_data2 <- read.delim(paste0("FUSION_final_SAVE_HSQ/",pop2,"_R2_hsq_4.tsv"))
    
    row_names <- sapply(strsplit(population_data1$GENE, "/"), `[`, 1)
    rownames(population_data1) <- row_names
    row_names <- sapply(strsplit(population_data2$GENE, "/"), `[`, 1)
    rownames(population_data2) <- row_names

    population_data1 <- population_data1[population_data1$HSQ.PV <= .05 & population_data1$HSQ < 1,]
    population_data2 <- population_data2[population_data2$HSQ.PV <= .05 & population_data2$HSQ < 1,]
    
    common_sites <- intersect(rownames(population_data1), rownames(population_data2))
    print(length(common_sites))
    population_data1 <- population_data1[rownames(population_data1) %in% common_sites,]
    population_data2 <- population_data2[rownames(population_data2) %in% common_sites,]
    output_file <- file(paste0(pop1, "_", pop2, "_fusion_snp_effects_final_save_hsq.txt"), "w")
    cat("Site\tModel\tSNP\teffect_pop1\teffect_pop2\n", file = output_file)
         for (sig_sites in rownames(population_data1)) {
                pop1_model <- population_data1[rownames(population_data1) == sig_sites, "MODEL"]
                pop2_model <- population_data2[rownames(population_data2) == sig_sites, "MODEL"]
                if(pop1_model == pop2_model)
                    {
                     print(sig_sites)
                     print(pop1_model)
                     if(pop1_model == "susie") {
                     load(paste0("FUSION_final_SAVE_HSQ/", pop1, "/", sig_sites, "/", sig_sites, ".wgt.RDat"))
                     pop1_wgt <- wgt.matrix
                     pop1_susie <- susiefit
                     load(paste0("FUSION_final_SAVE_HSQ/", pop2, "/", sig_sites, "/", sig_sites, ".wgt.RDat"))
                     pop2_wgt <- wgt.matrix
                     pop2_susie <- susiefit}
                    else{
                     load(paste0("FUSION_final_SAVE_HSQ/", pop1, "/", sig_sites, "/", sig_sites, ".wgt.RDat"))
                     pop1_wgt <- wgt.matrix
                     pop1_susie <- susiefit
                     load(paste0("FUSION_final_SAVE_HSQ/", pop2, "/", sig_sites, "/", sig_sites, ".wgt.RDat"))
                     pop2_wgt <- wgt.matrix
                     pop2_susie <- susiefit
                        }
                     best_model_wgts_pop1 <-  pop1_wgt[, colnames( pop1_wgt) == pop1_model]  
                     best_model_wgts_pop2 <-  pop2_wgt[, colnames( pop2_wgt) == pop2_model] 

                     if(pop1_model %in% c("enet", "lasso", "top1")){
                          best_model_wgts_pop1 <- best_model_wgts_pop1[abs(best_model_wgts_pop1)!= 0]
                          best_model_wgts_pop2 <- best_model_wgts_pop2[abs(best_model_wgts_pop2)!= 0]
                          print(length(best_model_wgts_pop1))
                          print(length(best_model_wgts_pop2))
                          overlapping_snps <- intersect(names(best_model_wgts_pop1), names(best_model_wgts_pop2))
                          print(length(overlapping_snps))
                          best_model_wgts_pop1 <- best_model_wgts_pop1[names(best_model_wgts_pop1) %in%  overlapping_snps]
                          best_model_wgts_pop2 <- best_model_wgts_pop2[names(best_model_wgts_pop2) %in%  overlapping_snps]
                          print(best_model_wgts_pop1)
                          print(best_model_wgts_pop2)
                          for(snp in overlapping_snps){
                             effect1 <- best_model_wgts_pop1[names(best_model_wgts_pop1) == snp][[1]]
                             effect2 <- best_model_wgts_pop2[names(best_model_wgts_pop2) == snp][[1]]
                             cat(sig_sites, pop1_model, snp, effect1, effect2, sep = "\t", file = output_file)
                             cat("\n", file = output_file)}
                          
                         }
                    else
                        {
                          best_model_wgts_pop1 <- pop1_susie$pip[pop1_susie$pip>.9]
                          best_model_wgts_pop2 <- pop2_susie$pip[pop2_susie$pip>.9]
                          print(length(best_model_wgts_pop1))
                          print(length(best_model_wgts_pop2))
                          overlapping_snps <- intersect(names(best_model_wgts_pop1), names(best_model_wgts_pop2))
                          print(length(overlapping_snps))
                          best_model_wgts_pop1 <- best_model_wgts_pop1[names(best_model_wgts_pop1) %in%  overlapping_snps]
                          best_model_wgts_pop2 <- best_model_wgts_pop2[names(best_model_wgts_pop2) %in%  overlapping_snps]
                          print(length(best_model_wgts_pop1))
                          print(length(best_model_wgts_pop2))
                          for(snp in overlapping_snps){
                             effect1 <- best_model_wgts_pop1[names(best_model_wgts_pop1) == snp][[1]]
                             effect2 <- best_model_wgts_pop2[names(best_model_wgts_pop2) == snp][[1]]
                             cat(sig_sites, pop1_model, snp, effect1, effect2, sep = "\t", file = output_file)
                             cat("\n", file = output_file)}
                              
                        }

                

                    print(head(best_model_wgts_pop1))
                    print(head(best_model_wgts_pop2))} }

             

 close(output_file)


