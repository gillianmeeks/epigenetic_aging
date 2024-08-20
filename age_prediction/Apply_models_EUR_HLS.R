
directory <- "LOOCV_test_train_split_final/"
files_meta_heritable <- list.files(directory, pattern = "*all_num_inds_138_age_transform=TRUEmeta_results_only=TRUEheritable_only=TRUE*", full.names = TRUE)
files_meta_not_heritable <- list.files(directory, pattern = "*all_num_inds_138_age_transform=TRUEmeta_results_only=TRUEheritable_only=FALSE*", full.names = TRUE)
files_not_meta_heritable <- list.files(directory, pattern = "*all_num_inds_138_age_transform=TRUEmeta_results_only=FALSEheritable_only=TRUE*", full.names = TRUE)
files_not_meta_not_heritable <- list.files(directory, pattern = "*all_num_inds_138_age_transform=TRUEmeta_results_only=FALSEheritable_only=FALSE*", full.names = TRUE)

load("final_age_predictor.RData")
for(df in c("meta_heritable_df", "meta_not_heritable_df", "not_meta_heritable_df", "not_meta_not_heritable_df"))
    {
    results <- get(df)
    print(df)
    print(paste("mean abs train error: ", round(mean(results$mean_abs_train_error),2)))
    #print(paste("mean real train error: ", round(mean(results$mean_real_train_error),2)))
    print(paste("mean abs test error: ", round(mean(results$mean_abs_test_african_error),2)))
    #print(paste("mean real test error: ", round(mean(results$mean_real_test_african_error),2)))   
    }
summary(meta_heritable_df$mean_abs_test_african_error)
median_meta_heritable_model <- meta_heritable_df[meta_heritable_df$mean_abs_test_african_error == 4.42,][1,1]


summary(meta_not_heritable_df$mean_abs_test_african_error)
median_meta_not_heritable_model <- meta_not_heritable_df[meta_not_heritable_df$mean_abs_test_african_error >= 4.05 & meta_not_heritable_df$mean_abs_test_african_error <= 4.06,][1,1]


summary(not_meta_not_heritable_df$mean_abs_test_african_error)
median_not_meta_not_heritable_model <- not_meta_not_heritable_df[not_meta_not_heritable_df$mean_abs_test_african_error >= 4.35 & not_meta_not_heritable_df$mean_abs_test_african_error <= 4.38,][1,1]


summary(not_meta_heritable_df$mean_abs_test_african_error)
median_not_meta_heritable_model <- not_meta_heritable_df[not_meta_heritable_df$mean_abs_test_african_error >= 4.82 & not_meta_heritable_df$mean_abs_test_african_error <= 4.84,][1,1]



undo_age_transformation <- function(transformed_age, adult_age) {
  original_age <- ifelse(transformed_age <= 0,
                          exp(transformed_age + log(adult_age + 1)) - 1,
                          transformed_age * (adult_age + 1) + adult_age)
  return(original_age)
}
adult_age <- 20
load("/share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/EMMAX/Horvath_test/datMethUsedNormalized_GSE78874.Rdata")
datSample$Ethnicity <- as.character(datSample$Ethnicity)
table(datSample$Ethnicity)
datSample <- datSample[, !colnames(datSample) %in% c("predictedAge", "age_pred_diff")]

online_results <- read.delim("noQC_no_pop_specific_missingness.csv", sep=",")
table(online_results$predictedTissue)
online_results <- online_results[online_results$Pop %in% c("Baka", "KHS", "European", "HispanicLatino", "Himba"),]
only_saliva <- online_results[online_results$predictedTissue %in% c("Saliva", "Blood PBMC", "Blood WB"), "SID"]

library(data.table)
library(glmnet)
library(ggplot2)

# Function to run the analysis
run_analysis <- function(files, pop, saliva, name) {
    abs_error_results_list <- list()
    real_error_results_list <- list()
    non_zero_list <- list()
    coefficients_list <- list()
    model_list <- list()
    cpg_list <- list()
  
  i <- 1
  
  if(saliva) {
    datSample <- datSample[rownames(datSample) %in% only_saliva, ]
  }
  print(table(datSample$Ethnicity))
  
  load("/share/hennlab/users/glmeeks/age_methylation/methylation_norm_EWAS/European_Hispanic_full_BMIQ.RData")
  
  if (any(is.na(European_Hispanic_BMIQ))) {
    print("There are NAs in the dataframe.")
  } else {
    print("There are no NAs in the dataframe.")
  }
  
  Caucasian_IDs <- rownames(datSample[datSample$Ethnicity == " Caucasian", ])
  HL_IDs <- rownames(datSample[datSample$Ethnicity == " Hispanic", ])
  HL_dat <- datSample[datSample$Ethnicity == " Hispanic", ]
  Caucasian_dat <- datSample[datSample$Ethnicity == " Caucasian", ]
  
  European_Hispanic_BMIQ <- na.omit(European_Hispanic_BMIQ)
  meth <- as.matrix(t(European_Hispanic_BMIQ))
  
  while (i <= length(files) && i < 101) {
    load(files[[i]])
    meth_ <- meth[, colnames(meth) %in% input_colnames]
    
    if (pop == "HSL") {
      meth_ <- meth_[rownames(meth_) %in% HL_IDs, ]
    } else {
      meth_ <- meth_[rownames(meth_) %in% Caucasian_IDs, ]
    }

    meth_ <- meth[, colnames(meth) %in% input_colnames]
    meth_ <- meth_[, c(input_colnames)]
    meth_ <- as.matrix(meth_)
    all_coefficients <- coeffs@Dimnames[[1]]

    ###order of columns needs to be the same as in the input model!!!!!!!!!!!!!!!!!!!
    predictions <- predict(best_model, s = models[[min_mse]]$lambda.min, newx = meth_)    
    predicted_undo <- undo_age_transformation(predictions, adult_age)
    
    if (pop == "HSL") {
      predicted_undo <- as.data.frame(predicted_undo[rownames(HL_dat), ])
    } else {
      predicted_undo <- as.data.frame(predicted_undo[rownames(Caucasian_dat), ])
    }
    
    colnames(predicted_undo) <- "predicted_age"
  
    if (pop == "HSL") {
      plot_data <- cbind(predicted_undo, HL_dat)
    } else {
      plot_data <- cbind(predicted_undo, Caucasian_dat)
    }
    plot_data$dif <- abs(plot_data$predicted_age - plot_data$Age)
    plot_data$real_dif <- plot_data$predicted_age - plot_data$Age 
    if("Inf" %in% plot_data$predicted_age ) {plot_data <- plot_data[- grep("Inf", plot_data$dif),]}else{plot_data<-plot_data}
    if("NaN" %in% plot_data$predicted_age ) {plot_data <- plot_data[- grep("NaN", plot_data$dif),]}else{plot_data<-plot_data}

    if("Inf" %in% plot_data$dif) {plot_data <- plot_data[- grep("Inf", plot_data$dif),]}else{plot_data<-plot_data}
    if("NaN" %in% plot_data$dif) {plot_data <- plot_data[- grep("NaN", plot_data$dif),]}else{plot_data<-plot_data}
  
    plot_data$dif <- abs(plot_data$predicted_age - plot_data$Age)
    plot_data$real_dif <- plot_data$predicted_age - plot_data$Age
    
    if ("Inf" %in% plot_data$predicted_age) {
      plot_data <- plot_data[-grep("Inf", plot_data$dif), ]
    }
    
    if ("Inf" %in% plot_data$dif) {
      plot_data <- plot_data[-grep("Inf", plot_data$dif), ]
    }
    
    mean_abs_error <- round(mean(plot_data$dif), 2)
    mean_real_error <- round(mean(plot_data$real_dif), 2)
    
    abs_error_results_list[[i]] <- mean_abs_error
    real_error_results_list[[i]] <- mean_real_error
    non_zero_list[[i]] <- non_zero
    coefficients_list[[i]] <- coeffs
    model_list[[i]] <- files[[i]]
    cpg_list[[i]] <- length(non_zero) - 1
    i <- i + 1
  }
  
  results_df <- data.table(model = unlist(model_list), num_cpgs = unlist(cpg_list), mean_abs_error = unlist(abs_error_results_list), mean_real_error = unlist(real_error_results_list))
  results_df <- as.data.frame(results_df)
  
  if (saliva) {
    assign(paste0(name, "_applied_to_", pop, "_saliva"), results_df, envir = .GlobalEnv)
  } else {
    assign(paste0(name, "_applied_to_", pop, "_allsamps"), results_df, envir = .GlobalEnv)
  }
}

# Parameters
file_lists <- list("files_not_meta_not_heritable", "files_meta_not_heritable", "files_not_meta_heritable", "files_meta_heritable")
pops <- c("Caucasian", "HSL")
saliva_vals <- c(TRUE, FALSE)

# Main loop
for (name in file_lists) {
  files <- get(name)
  
  for (pop in pops) {
    for (saliva in saliva_vals) {
      run_analysis(files, pop, saliva, name)
    }
  }
}

save(list=ls(), file="models_to_HSL_Euro.RData")
