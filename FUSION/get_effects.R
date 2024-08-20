args <- commandArgs(trailingOnly = TRUE)
pop <- as.character(args[1])
#pop <- "KHS"

#need to re-run it to not just include snps with the same model in both populations.
    population_data <- read.delim(paste0("FUSION_results_FINAL/",pop,"_R2_hsq_4.tsv"))
    head(population_data)
     row_names <- sapply(strsplit(population_data$GENE, "/"), `[`, 1)
    rownames(population_data) <- row_names  
    population_data <- population_data[population_data$HSQ.PV <= .05 & population_data$HSQ < 1 & population_data$HSQ > 0,]

    output_file <- file(paste0(pop, "_fusion_snp_effects_final_save_hsq.txt"), "w")
    cat("Site\tModel\t\tSNP\teffect\n", file = output_file)
         for (sig_sites in rownames(population_data)) {
             
                pop_model <- population_data[rownames(population_data) == sig_sites, "MODEL"]
                print(pop_model)
                    load(paste0("FUSION_final_SAVE_HSQ/", pop, "/", sig_sites, "/", sig_sites, ".wgt.RDat"))
                     pop_wgt <- wgt.matrix
                     if(pop_model == "susie") { 
                      pop_susie <- susiefit
                      best_model_wgts_pop <- pop_susie$pip[pop_susie$pip>.9]
                     print(best_model_wgts_pop )
                     }
                     if(pop_model %in% c("enet", "lasso", "top1")){
                        best_model_wgts_pop <-  pop_wgt[, colnames( pop_wgt) == pop_model]  
                        best_model_wgts_pop <- best_model_wgts_pop[abs(best_model_wgts_pop)!= 0]
                         print( best_model_wgts_pop )
                         }
                         # Iterate over the named effects and write each one as a row
                      for (snp in names(best_model_wgts_pop)) {
                        effect <- best_model_wgts_pop[snp]
                        cat(sig_sites, pop_model, snp, effect, sep = "\t", file = output_file, append = TRUE)
                        cat("\n", file = output_file, append = TRUE)
                      }
                    }
         

 close(output_file)

