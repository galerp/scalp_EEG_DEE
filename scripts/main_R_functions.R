library(ggpubr)
library(tidyverse)
library(reshape2)
library(reticulate)
library(Hmisc)
library(caret)
library(randomForest)
library(pROC)

##### Age matching Function #####


# This function performs age matching from a disease cohort "df_gene" and a 
# control cohort "df_control". Note: it assumes "df_gene" has more datapoints
# Each control data point can only be used once.

# Input:
#      df_gene: A dataframe from a cohort with a specific genetic disorder 
#               in a wide format containing EEG features for each data point.
#               Each row is one EEG from one individual.
#      df_control: A dataframe from a large control cohort in a wide format 
#               containing EEG features for each data point.
#               Each row is one EEG from one individual.
# Output:
#      matched_control_rows: A subset of "df_control" that contains the 
#               age matched controls for the "df_gene" cohort.

find_age_matched_controls <- function(df_gene, df_control) {
  
  # Add an index column to each dataframe for tracking
  df_gene <- df_gene %>% mutate(index_gene = row_number())
  df_control <- df_control %>% mutate(index_control = row_number())
  
  # Initialize a list to hold the matched controls
  matched_controls <- vector("list", nrow(df_gene))
  
  # Copy of the control dataframe to modify during matching
  available_controls <- df_control
  
  # Loop through each individual in the gene group to find the closest age match from available controls
  for (i in seq_along(matched_controls)) {
    current_age <- df_gene$age[i]
    
    # Find the control with the closest age from the available pool
    closest <- available_controls %>%
      mutate(age_diff = abs(age - current_age)) %>%
      arrange(age_diff) %>%
      slice_head(n = 1)
    
    # Store the matched control index
    matched_controls[[i]] <- closest$index_control
    
    # Remove the matched control from the available pool
    available_controls <- available_controls %>%
      filter(index_control != closest$index_control)
  }
  
  # Use the matched control indices to extract the corresponding rows from df_control
  matched_control_rows <- df_control %>%
    filter(index_control %in% unlist(matched_controls)) %>% 
    select(-index_control)
  
  return(matched_control_rows)
}


##### RANDOM FOREST #####

# This function trains random forest models for every gene in a dataframe.
# It first uses the find_age_matched_controls to perform age matching for 
# every individual in the gene group, gathering from the larger Control group
# It performs leave one out cross validation (LOOCV).

# Input:
#      data_pretrain: A dataframe in a wide format containing features for each
#                     individual. The first four columns are patient identifiers,
#                     age, and gene group. If there are columns after that,
#                     they contain spectral information. Each row is one EEG
#                     from one individual.
# Output:
#      results_rf: A dataframe containing the performance and feature importance
#                  from every iteration of the LOOCV from each gene.

# This handles if it is taking in data that just has age as a feature

rf_gene_control <- function(data_pretrain) {
  
  if(length(data_pretrain)>4){
    allfeats = data_pretrain[,5:length(data_pretrain)]
    nfeats = length(allfeats)
    feat_names = names(allfeats)
    
    results_df <- as.data.frame(matrix(nrow=0, ncol=nfeats+7))
    names(results_df) <- c(c("gene","Iteration", "prediction_prob", "prediction", 
                             "age_case","actual_outcome", "age"), feat_names)
  }else{
    results_df <- as.data.frame(matrix(nrow=0, ncol=7))
    names(results_df) <- c("gene","Iteration", "prediction_prob", "prediction", 
                           "age_case","actual_outcome", "age")
    
  }
  
  
  genes = unique(data_pretrain$gene)
  genes = genes[genes != "Controls"]
  
  for (g in genes) {
    data_gene <- data_pretrain %>% 
      filter(gene == g)
    
    df_control <- data_pretrain %>%
      filter(gene == "Controls")
    
    data_combined <- find_age_matched_controls(data_gene, df_control)
    data_combined <- rbind(data_gene, data_combined) %>%
      mutate(gene_binary = ifelse(gene == g, 1, 0))
    data_combined$gene_binary <- factor(data_combined$gene_binary)
    
    # Perform LOOCV
    for (i in 1:nrow(data_combined)) {
      training <- data_combined[-i, ]
      testing <- data_combined[i, ]
      
      # Train Random Forest model
      rf_model <- randomForest(gene_binary ~ ., data = training %>% select(-c(patient, gene, age_group)), ntree = 500, importance = TRUE)
      
      # Extract feature importance
      importance_vals <- importance(rf_model) %>% as.data.frame()
      
      # Predict and evaluate
      prediction_prob <- predict(rf_model, newdata = testing %>% select(-c(patient, gene, age_group)), type = "prob")[,2]  # Get probability of class 1
      prediction <- as.numeric(prediction_prob > 0.5)
      actual <- testing$gene_binary
      
      # Prepare the results dataframe
      # For MeanDecreaseGini
      cur_df <-  as.data.frame(t(c(g, i, testing$age[1],prediction_prob, prediction, actual, importance_vals$MeanDecreaseGini)))
      # for MeanDecreaseAccuracy
      # cur_df <-  as.data.frame(t(c(g, i, testing$age[1],prediction_prob, prediction, actual, importance_vals$MeanDecreaseAccuracy)))
      
      names(cur_df) <- c("gene", "Iteration", "age_case", "prediction_prob", "prediction", "actual_outcome", 
                         rownames(importance_vals))
      
      # Store results
      results_df <- rbind(results_df, cur_df)
    }
  }
  results_rf = results_df %>% 
    mutate(actual_outcome = as.numeric(actual_outcome)-1) %>%
    mutate(prediction_prob = as.numeric(prediction_prob)) %>% 
    mutate(hit = case_when(prediction == actual_outcome ~ 1,
                           TRUE ~ 0))%>%
    mutate(FP = case_when(prediction == 1 & actual_outcome == 0 ~ 1,
                          TRUE ~ 0)) %>%
    mutate(TP = case_when(prediction == 1 & actual_outcome == 1 ~ 1,
                          TRUE ~ 0)) %>%
    mutate(TN = case_when(prediction == 0 & actual_outcome == 0 ~ 1,
                          TRUE ~ 0))%>%
    mutate(FN = case_when(prediction == 0 & actual_outcome == 1 ~ 1,
                          TRUE ~ 0))
  return(results_rf)
}



##### Seizure frequence Random forest with LOOCV #####

# Assumes balanced dataset
seizfreq_rf_loocv <- function(train_df) {
  
  results_df <- data.frame()
  
  for (i in 1:nrow(train_df)) {
    
    # print(i)
    training <- train_df[-i, ]
    testing <- train_df[i, ]# Separate training and testing data
    
    
    # Train model
    rf_model <- randomForest(seiz_binary ~ ., data = training %>% select(-c(patient,seiz_freqs,seiz_age)), ntree = 500, importance = TRUE)
    prediction_prob <- predict(rf_model, newdata = testing %>% select(-c(patient, seiz_freqs,seiz_age)), type = "prob")[,2]
    prediction <- as.numeric(prediction_prob > 0.5)
    
    # Extract feature importance
    importance_vals <- importance(rf_model) %>% as.data.frame()
    feat_import <- as.data.frame(t(importance_vals$MeanDecreaseGini))
    names(feat_import) = rownames(importance_vals)
    
    cur_df <-  cbind(testing$age,testing$patient,prediction_prob, prediction, as.numeric(testing$seiz_binary[1])) %>% as.data.frame()
    names(cur_df) <- c("age_eeg", "patient","prediction_prob", "prediction", "actual_outcome")
    cur_df$Iteration <- i
    
    #####  Null models ----
    null_model <- randomForest(seiz_binary ~ age, data = training %>% select(age, seiz_binary), ntree = 500, importance = TRUE)
    null_predictions_prob <- predict(null_model, newdata = testing %>% select(age), type = "prob")[,2]
    null_predictions <- as.numeric(null_predictions_prob > 0.5)
    
    
    ##### ----
    null_cur_df <- cbind(null_predictions_prob, null_predictions) %>% as.data.frame()
    
    names(null_cur_df) <- c("null_prediction_prob", "null_prediction")
    # null_cur_df$Iteration <- i
    cur_df <- cur_df %>% cbind(null_cur_df)%>% cbind(feat_import)
    
    results_df <- rbind(results_df, cur_df)
    # results_df <- rbind(results_df, null_cur_df)
  }
  
  return(results_df)
  
}



##### GMFM Random forest with bootstraping #####

# This function trains random forest models predicting GMFM scores.
# It runs 1000 bootstraps wiht an 80/20 training-testing split. It trains and
# tests an EEG model (with EEG features) and a null model with just age as a
# feature.

# Input:
#      train_df: A dataframe in a wide format containing features for each
#                     individual. The first four columns are patient identifiers,
#                     age, and gene group. If there are columns after that,
#                     they contain spectral information. Each row is one EEG
#                     from one individual.
# Output:
#      results_rf: A dataframe containing the performance and feature importance
#                  from every iteration of the bootstrap. Performance is 
#                  is primarily evaluated with RMSE and RMSE.

# This handles if it is taking in data that just has age as a feature

gmfm_rf <- function(train_df) {
  
  
  results_df <- data.frame()
  null_results_df <- data.frame()
  n_train <- round(0.8*nrow(train_df))
  
  
  for (i in 1:1000) {
    
    train_rows <- sample(seq(1,nrow(train_df),1),n_train)
    training <- train_df[train_rows,]
    testing <- train_df[-train_rows,]  # Separate training and testing data
    
    
    # Train model
    rf_model <- randomForest(MEAS_VALUE ~ ., data = training %>% select(-c(patient)), ntree = 500, importance = TRUE)
    predictions <- predict(rf_model, newdata = testing %>% select(-c(patient)), type = "response")
    # Extract feature importance
    importance_vals <- importance(rf_model) %>% as.data.frame()
    feat_import <- as.data.frame(t(importance_vals$IncNodePurity))
    names(feat_import) = rownames(importance_vals)

    #### Null models ----

    # linear_model <- lm(MEAS_VALUE ~ age_test, data = training)
    # null_predictions <- predict(linear_model, newdata = testing %>% select(age_test), type = "response")
    null_model <- randomForest(MEAS_VALUE ~ ., data = training %>% select(age_test, MEAS_VALUE), ntree = 500, importance = TRUE)
    null_predictions <- predict(null_model, newdata = testing %>% select(age_test), type = "response")
    
    #### ----
    
    
    actual <- testing$MEAS_VALUE
    
    
    # Prepare the results dataframe
    cur_df <-  as.data.frame(matrix(c(predictions,null_predictions,actual),nrow = length(actual)))
    names(cur_df) <- c( "prediction", "null_predictions", "actual_outcome")
    cur_df$iteration = i
    cur_df$age = testing$age_test
    
    rmse = sqrt(mean((actual - predictions)^2))
    null_rmse = sqrt(mean((actual - null_predictions)^2))
    mae = mean(abs(actual - predictions))
    null_mae = mean(abs(actual - null_predictions))
    
    rmse_dif = null_rmse - rmse
    mae_dif = null_mae - mae
    
    cur_df$mae = mae
    cur_df$null_mae = null_mae
    cur_df$mae_dif = mae_dif
    
    cur_df$rmse = rmse
    cur_df$null_rmse = null_rmse
    cur_df$rmse_dif = rmse_dif
    
    cur_df$dif_comp = 0
    
    cur_df = cur_df %>% cbind(feat_import)
    # Store results
    results_df <- rbind(results_df, cur_df)
    # null_results_df <- rbind(null_results_df, null_cur_df)
    
    
  }
  return(results_df)
  
}
