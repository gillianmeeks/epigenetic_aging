
start_time <- Sys.time()

library(glmnet)
library(ggplot2)
library(foreach)

print("loaded")


  # Arguments
args <- commandArgs(trailingOnly = TRUE)
task <- as.numeric(args[1])
meta_results <- FALSE
heritable_only <- FALSE


age_transform <- TRUE
## saliva <- as.logical(args[3])

pop <- "all"


print("here")
      load("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/meqtl_cpgs.RData")
      load("/share/hennlab/users/glmeeks/age_methylation/EMMAX_meqtl/heritable_cpgs.RData")
      
      # Additional file loading checks
      required_files <- c(
        "/share/hennlab/users/glmeeks/age_methylation/age_prediction/merged_data_all_phenos.RData",
        "/share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/Diverse_Age_QC_methylation_phenos.RData",
        "/share/hennlab/users/glmeeks/age_methylation/age_prediction/all_all_sites_no_regressed.RData",
        "/share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/EMMAX/metagen_all3_no_meqtl_regress.RData",
        "/share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/European_Hispanic_full_BMIQ.RData"
      )
      for (file in required_files) {
        if (!file.exists(file)) {
          stop(paste("File not found:", file))
        }
        load(file)
      }
    meth <- as.data.frame(no_residuals_df)
    
      ######### Data Preparation #########
    set.seed(NULL)
    Baka_inds <- c("B01E", "B03F", "B03H", "B04E", "B05E", "B06E", "B08E", "B09F", "B09H", "B12E", "B13F", "B13H",
                       "B16E", "B14F", "B14H", "B17E", "B19F", "B19H", "B20E","B22F", "B22H", "B23F", "B23H", 
                   "B03E", "B05F", "B05H", "B09E", "B12F", "B12H", "B13E", "B16H", "B19E", "B20F", "B20H", "B22E")
    #randomly shuffle so the split is different each time
    Baka_inds <- sample(Baka_inds, length(Baka_inds))
    #get prefix and child (E) or parent (H or F)
    first_three <- substr(Baka_inds, 1, 3)
    last_char <- substr(Baka_inds, 4, 4)
    df <- data.frame(ID = Baka_inds, First_Three = first_three, Last_Char = last_char)
    df <- df[order(match(df$ID, Baka_inds)), ]
    #if duplicated first three then a member of a trio
    duplicates <- df[duplicated(df$First_Three), "First_Three"]
    df$trio <- NA
    df[df$First_Three %in% duplicates, "trio"] <- TRUE
    df[!df$First_Three %in% duplicates, "trio"] <- FALSE
    #two that are just parents no egos so just put in no trios
    df[df$First_Three %in% c("B23", "B14"), "trio"]  <- FALSE
    #sample 7 no trios names
    no_trio_for_train <- sample(df[!df$trio, "ID"], 7)
    #take remainng no trios for test
    no_trio_for_test <- df[!df$trio & !df$ID %in% no_trio_for_train, "ID"]
    #B16 is just one parent one child and just always put the parent in one and child in the other (based on the j value)
    trios <- df[df$trio & !df$First_Three == "B16" ,]
    trio_prefixes <- sample(unique(trios$First_Three), length(unique(trios$First_Three)))
    trio_prefixes 
    parents_train <- trio_prefixes[c(1,2,3,4,5,6)]
    parents_test <- trio_prefixes[!trio_prefixes %in% parents_train]
    
    parentstrainF <- c()
    parentstrainH <- c()
    for (prefix in parents_train) {
      parentstrainF <- c(parentstrainF, paste0(prefix, "F"))
      parentstrainH <- c(parentstrainH, paste0(prefix, "H"))
    }
      parentstrain <- c(parentstrainF, parentstrainH)
    
    parentstestF <- c()
    parentstestH <- c()
    for (prefix in parents_test) {
      parentstestF <- c(parentstestF, paste0(prefix, "F"))
      parentstestH <- c(parentstestH, paste0(prefix, "H"))
    }
    parentstest <- c(parentstestF, parentstestH)
    
    kids_train <- trio_prefixes[c(7,8)]
    kids_test <- trio_prefixes[c(1,2,3,4,5,6)]
    
    kidstrainE <- c()
    for (prefix in kids_train) {
      kidstrainE <- c(kidstrainE, paste0(prefix, "E"))
    }
    kidstrain <- c(kidstrainE)
    
    kidstestE <- c()
    for (prefix in kids_test) {
      kidstestE <- c(kidstestE , paste0(prefix, "E"))
    }
      kidstest <- c(kidstestE)
    
    #put the 1parent child only separate
    j <- runif(1)
    if(j %% 2 == 0) {training <- c(kidstrain, parentstrain, "B16E", no_trio_for_train)    
                     testing <- c(kidstest, parentstest, "B16H", no_trio_for_test)
    } else {training <- c(kidstrain, parentstrain, "B16H", no_trio_for_train)    
                          testing <- c(kidstest, parentstest, "B16E", no_trio_for_test)}
    print(paste0("length of Baka test: ", length(testing), "length of Baka train: ", length(training)))
    ####Himba and KHS splits to 70% train
    Himba_train_split <- round(51*.7,0)
    KHS_train_split <- round(52*.7,0)
    
    Himba <- colnames(Himba_merged)
    print(length(Himba[Himba %in% rownames(meth)]))
    print(length(Himba[Himba %in% rownames(Himba_pheno_merged)]))
        #shuffle
    Himba <- sample(Himba, length(Himba))
    Himba <- Himba[!Himba %in% c("HMB181_2")]
    KHS <- colnames(KHS_merged)
        #shuffle
    KHS <- sample(KHS, length(KHS))
    
    Himba_train<- sample(Himba, Himba_train_split)
    KHS_train<- sample(KHS, KHS_train_split)
    Himba_test <- Himba[!Himba %in% Himba_train]
    KHS_test <- KHS[!KHS %in% KHS_train]
    training <- c(training, Himba_train, KHS_train)
    testing <- c(testing, Himba_test, KHS_test)
    
    undo_age_transformation <- function(transformed_age, adult_age) {
      original_age <- ifelse(transformed_age <= 0,
                              exp(transformed_age + log(adult_age + 1)) - 1,
                              transformed_age * (adult_age + 1) + adult_age)
      return(original_age)
    }

      age <- pheno_merged[,c("ID", "age")]
      age.train <- pheno_merged[pheno_merged$ID %in% training,]
      age.test <- pheno_merged[pheno_merged$ID %in% testing,]
      meth <- meth[, apply(meth, 2, function(x) all(!is.na(x)))]
      online_results <- read.delim("/share/hennlab/users/glmeeks/age_methylation/age_prediction/noQC_no_pop_specific_missingness.csv", sep = ",")
      online_results <- online_results[online_results$Pop %in% c("Baka", "KHS", "European", "HispanicLatino", "Himba"),]
      online_results[online_results$SID == "HMB494.2", "SID"] <- "HMB494-2"
                           ##shouldn't do only saliva because then I have to hand split the 4 Baka that aren't predicted saliva
      # only_saliva <- online_results[online_results$predictedTissue == "Saliva", "SID"]
      # if (saliva) {
      #   meth <- meth[rownames(meth) %in% only_saliva,]
      #   pheno_merged <- pheno_merged[pheno_merged$ID %in% only_saliva,]
      # } else {
      #   meth <- meth
      # }
                           #same cpgs as Euro
      European_Hispanic_BMIQ <- na.omit(European_Hispanic_BMIQ)                     
      meth <- meth[, colnames(meth) %in% rownames(European_Hispanic_BMIQ)]

      sig <- subset(meta_gen, P.value.fix < .05 / nrow(meta_gen))
     
      any_heritable_or_meqtl <- unique(c(pop_spec_meqtl_cpgs, pop_spec_heritable_cpgs))
      #any_heritable_or_meqtl <- unique(pop_spec_meqtl_cpgs)
      if (heritable_only) {
        meth <- meth[, colnames(meth) %in% any_heritable_or_meqtl]
      } else {
        meth <- meth[, !colnames(meth) %in% any_heritable_or_meqtl]
      }

      if (meta_results) {
        meth <- meth[, colnames(meth) %in% sig$CPG.Labels]
      } else {
        meth <- meth
      }

      input_colnames <- colnames(meth)
       ####test subset####                    
     ####test small subset 
    #meth<- meth[,1:1000]    
     ############                    
    ncol(meth)  
    meth <- meth[c(rownames(pheno_merged)),]
    training <- meth[rownames(meth) %in% training,]
    testing <- meth[rownames(meth) %in% testing,]
    ########age transform####
    adult_age <- 20
      
    # if(log_){meth <- log(meth)}
    # else{meth <- meth}
    print(head(meth[,1:6]))
                             
    train <- as.matrix(training, )
    test <- as.matrix(testing, )
    
    age.train_ <- age.train[,"age"]
    age.test_ <- age.test[,"age"]
    
    
    if (age_transform) {
      # Apply the Horvath age transform logic to age.train
      age.train <- ifelse(age.train_ <= adult_age,
                          log(age.train_ + 1) - log(adult_age + 1),
                          (age.train_ - adult_age) / (adult_age + 1))
      
      age.test <- ifelse(age.test_ <= adult_age,
                         log(age.test_ + 1) - log(adult_age + 1),
                         (age.test_ - adult_age) / (adult_age + 1))
    } else {
      age.train <- age.train_
      age.test <- age.test_
    }
    models <- list()

    for (i in 0:20) {
      name <- paste0("alpha", i/20)
      models[[name]] <- cv.glmnet(train, age.train, type.measure = "mse", grouped=FALSE, keep=TRUE, nfolds=nrow(train), alpha = i/20, family = "gaussian")
     
    }

    results <- data.frame()
                     ###best model based on lowest mean squared error in test data
 for (i in 0:20) {
      name <- paste0("alpha", i/20)
      ## Use each model to predict 'y' given the Testing dataset 
      predicted <- predict(models[[name]], s = models[[name]]$lambda.min, newx = test)
      ## Ensure dimensions match
      age.test <- as.vector(age.test)
      predicted <- as.vector(predicted)

      ## Calculate the Mean Squared Error...
      mse <- mean((age.test - predicted)^2)

      ## Store the results
      temp <- data.frame(alpha = i/20, mse = mse, name = name)
      results <- rbind(results, temp)
    }

# Print the alpha value and corresponding MSE of the best model (best alpha model, and then print the lambda that had the lowest mean se across all the folds)

min_mse <- results[results$mse == min(results$mse), "name"][1]  
cat("Best Alpha:", min_mse, "\n")# Select the first model with the minimum MSE
coeffs <- predict(models[[min_mse]], type = "coef")
non_zero_coefficients <- coeffs[coeffs[, 1] != 0, ]
print(paste0("number of nonzero coeffs: ", length(non_zero_coefficients)))
non_zero <- names(non_zero_coefficients)
       
   
    if(age_transform){
        # Undo the Horvath age transform on training, test, and predicted ages
        age_train_undo <- undo_age_transformation(age.train, adult_age)
        age_test_undo <- undo_age_transformation(age.test, adult_age)
        # Predict the training data
        predicted_train <- predict(models[[min_mse]], s = models[[min_mse]]$lambda.min, newx = train)
        # Undo the horvath age transform on training predicted ages
        predicted_train_undo <- undo_age_transformation(predicted_train, adult_age)
        predicted_best <- predict(models[[min_mse]], s = models[[min_mse]]$lambda.min, newx = test)
        # Store the predictions in the predicted_best variable
        predicted_best <- as.vector(predicted_best)
        # Undo the Horvath age transform on predicted ages
        predicted_undo <- undo_age_transformation(predicted_best, adult_age)

    } else{
         age_train_undo <- age.train
         age_test_undo <- age.test
         predicted_train <- predict(models[[min_mse]], s = models[[min_mse]]$lambda.min, newx = train)
         predicted_train_undo <- predicted_train
         predicted_best <- predict(models[[min_mse]], s = models[[min_mse]]$lambda.min, newx = test)
         predicted_best <- as.vector(predicted_best)
         predicted_undo <- predicted_best
    }
   
    train <- data.frame(
      RealAge = as.vector(age_train_undo),
      PredictedAge = as.vector(predicted_train_undo),
      Dataset = "Training"
    )

    test <- data.frame(
      RealAge = as.vector(age_test_undo),
      PredictedAge = as.vector(predicted_undo),
      Dataset = "Testing"
    )

    # Combine the data frames
    plot_data <- rbind(train, test)
    plot_data$dif <- abs(plot_data$RealAge - plot_data$PredictedAge)
    plot_data$real_dif <- plot_data$PredictedAge - plot_data$RealAge
    abs_error_train <- mean(plot_data[plot_data$Dataset == "Training", "dif"], na.rm = TRUE)
    abs_error_test <- mean(plot_data[plot_data$Dataset == "Testing", "dif"], na.rm = TRUE)
    real_error_train <- mean(plot_data[plot_data$Dataset == "Training", "real_dif"], na.rm = TRUE)
    real_error_test <- mean(plot_data[plot_data$Dataset == "Testing", "real_dif"], na.rm = TRUE)                 

    # Plot real vs predicted ages using ggplot2
k <- ggplot(plot_data, aes(x = RealAge, y = PredictedAge, color = Dataset)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    x = "Real Age", 
    y = "Predicted Age",
    title = paste(
      pop, 
      "Real vs Predicted transform =", age_transform, 
      "meta_results_only =", meta_results, "\n", 
      "heritable_only =", heritable_only, "\n",
      "number_cpgs =", ncol(meth), "\n",
      "Training abs Mean Error:", round(abs_error_train, 2), " years", "\n",
      "Testing abs Mean Error:", round(abs_error_test, 2), " years", "\n",
      "Correlation real and predicted:", cor(plot_data$RealAge, plot_data$PredictedAge)
    )
  ) +
  theme_minimal() +
  scale_color_discrete(name = "Dataset")  # Add a legend



file_name <- paste0("/share/hennlab/users/glmeeks/age_methylation/age_prediction/LOOCV_test_train_split_final_svgs/", pop, "_age_transform=", age_transform, "meta_results_only=", meta_results, "heritable_only=", heritable_only, "_num_inds_", nrow(meth),
                    "_number_cpgs=", ncol(meth), "_TrainingMeanError=", round(abs_error_train, 2),
                    "_Corr=", round(cor(plot_data$RealAge, plot_data$PredictedAge), 2), ".svg")

  # Calculate correlation between real and predicted ages
  correlation <- cor(plot_data$RealAge, plot_data$PredictedAge)
  
  # Create a data frame with the results for this iteration
  iteration_results <- data.frame(
    #Iteration = iteration,
    Abs_Error_Train = abs_error_train,
    Abs_Error_Test = abs_error_test,
    Real_Error_Train = real_error_train,
    Real_Error_Test = real_error_test,    
    Correlation = correlation
  )

    print(file_name)


                     # Save the SVG files
ggsave(file_name, plot = k, width = 10, height = 8)
                     lambda <- models[[min_mse]]$lambda.1se
                     best_model <- models[[min_mse]]
save(list=c("non_zero", "models", "min_mse", "best_model", "coeffs", "input_colnames", "plot_data"), file= paste0("/share/hennlab/users/glmeeks/age_methylation/age_prediction/LOOCV_test_train_split_final/", task, pop, "_num_inds_", nrow(meth), "_age_transform=", age_transform, "meta_results_only=", meta_results, "heritable_only=", heritable_only,
                    "_number_cpgs=", ncol(meth), "_TrainingMeanError=", round(abs_error_train, 2),
                    "_Corr=", round(cor(plot_data$RealAge, plot_data$PredictedAge), 2), ".RData"))


                      
# write.csv(iteration_results,file= paste0(task,"_cross_val_train_split_with_LOOCV_age_transform=", age_transform, "meta_results_only=", meta_results, "heritable_only=", heritable_only, "_num_inds_", nrow(meth), ".csv"), row.names = FALSE)


end_time <- Sys.time()                            
runtime <- end_time - start_time

runtime_hours <- as.numeric(runtime, units = "hours")
cat("Runtime in hours:", runtime_hours, "\n") 


                           
