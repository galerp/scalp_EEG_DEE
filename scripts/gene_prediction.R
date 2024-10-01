###################
# Description: 
# This file creates the random forest models to predict both genes from controls
# and genes from each other. These results are seen in Figure 4 of Galer et al.
# Note, results may differ very slightly due to the stochastic nature of random 
# forest models. For example the top 8 features may shift slightly in order 
# (ordered by median). Also, note LOOCV with close age matching and 1:1 ratio
# between both classes results in perfect or near perfect AUCs for models trained on 
# only age due to majority class prediction (AUC prediction probability is flipped).

###################

# Set working directory
setwd("/Volumes/helbig_lab/Users/galerp/EEG/manuscript/scalp_EEG_DEE/")

# Loads primary functions
source("scripts/main_R_functions.R")

# Spectral information for every electrode
psd_bp_rel <- read_csv("data/psd_bp_gene_controls.csv")

####################
# Median Across Scalp
####################

psd_bp_rel_pat <- psd_bp_rel %>% 
  group_by(patient, age, age_group, gene, bandpower) %>% 
  dplyr::summarise(med_feat = median(power)) %>% 
  filter(age_group !="Neonates")

####################
# Median Across Lobes
####################

psd_bp_rel_lobe <- psd_bp_rel %>% 
  group_by(patient, age, age_group, gene, lobe, bandpower) %>% 
  dplyr::summarise(med_feat = median(power)) %>% 
  filter(age_group !="Neonates")


##############################################
# Combine all Features
##############################################

# Just Spatial Features
full_feat_lobe <- psd_bp_rel_lobe %>%
  ungroup() %>% 
  mutate(feature = paste(bandpower, lobe, sep = "_")) %>% 
  mutate(feature = gsub("_power_rel","",feature)) %>% 
  mutate(feature = gsub("Frontal","Frnt",feature)) %>% 
  mutate(feature = gsub("Occipital","Occ",feature)) %>% 
  mutate(feature = gsub("Temporal","Tmp",feature)) %>% 
  mutate(feature = gsub("Central","Cnt",feature)) %>% 
  mutate(feature = gsub("Parietal","Par",feature)) %>% 
  select(-lobe) 



# Spatial and Non-Spatial Features
full_feat_all <- psd_bp_rel_pat %>% 
  ungroup() %>% 
  mutate(feature = gsub("_power_rel","",bandpower)) %>% 
  rbind(full_feat_lobe)

########
# Covert Features DF to Wide Format
########

feat_lobe_wide <- full_feat_lobe %>%
  pivot_wider(names_from = feature, values_from = med_feat,
              names_prefix = "bp_", id_cols = c(patient, gene, age, age_group)) 

# Change feature names for easier viewing
wide_names = names(feat_lobe_wide)
new_names = gsub("_power","",wide_names)
new_names = gsub("bp_","",new_names)
new_names = gsub("_all","",new_names)
names(feat_lobe_wide) = new_names

feat_all_wide <- full_feat_all %>%
  pivot_wider(names_from = feature, values_from = med_feat,
              names_prefix = "bp_", id_cols = c(patient, gene, age, age_group)) 
# Change feature names for easier viewing
wide_names = names(feat_all_wide)
new_names = gsub("_power","",wide_names)
new_names = gsub("bp_","",new_names)
new_names = gsub("_all","",new_names)
names(feat_all_wide) = new_names


feat_age_wide <- full_feat_all %>%
  select(patient, gene, age, age_group) %>% 
  distinct()

################################
# Train Models
################################

################
# Just Lobe Spectral Features
################


data_pretrain = feat_lobe_wide %>%
  ungroup() 

results_main_rf <- rf_gene_control(data_pretrain)

################
# Only Non-spatial Features
################
all_feat_names = names(feat_all_wide)
spatial_substrings = c("Cnt","Par","Occ","Tmp","Frnt")

# Find indices of all names that contain any of the substrings
matched_indices <- unique(unlist(sapply(spatial_substrings, function(sub) grep(sub, all_feat_names))))

# Exclude these indices from the original list
columns_to_keep <- all_feat_names[-matched_indices]

data_pretrain = feat_all_wide %>%
  ungroup() %>% 
  select(all_of(columns_to_keep))

results_nospace_rf <- rf_gene_control(data_pretrain)

################
# All Spectral Features
################

data_pretrain = feat_all_wide %>%
  ungroup() 

results_all_rf <- rf_gene_control(data_pretrain)

####################################
# STXBP1 Results
####################################

stx_results = results_main_rf %>% 
  filter(gene == "STXBP1") %>% 
  mutate(prediction_prob = as.numeric(prediction_prob))


stx_accuracy = sum(stx_results$hit)/nrow(stx_results)
stx_precision = sum(stx_results$TP)/(sum(stx_results$TP)+sum(stx_results$FP))
stx_recall = sum(stx_results$TP)/(sum(stx_results$TP)+sum(stx_results$FN))
stx_f1 = (2 * stx_precision * stx_recall)/(stx_precision+stx_recall)

print(paste0("Accuracy: ", stx_accuracy))

print(paste0("Precision: ", stx_precision))

print(paste0("Recall: ", stx_recall))

print(paste0("F1: ", stx_f1))


######
# Get ROC and AUC
######

# Main Results
# Results with no spatial features
stx_results_nospace <- results_nospace_rf %>% 
  filter(gene == "STXBP1") %>% 
  mutate(prediction_prob = as.numeric(prediction_prob))
roc_data_stx_rf_nospace<- roc(stx_results_nospace$actual_outcome, stx_results_nospace$prediction_prob)
auc_value <- auc(roc_data_stx_rf_nospace)
print(auc_value)


# Main results with localized features
stx_results = results_main_rf %>% 
  filter(gene == "STXBP1") %>% 
  mutate(prediction_prob = as.numeric(prediction_prob))
roc_data_stx_rf <- roc(stx_results$actual_outcome, stx_results$prediction_prob)
auc_value <- auc(roc_data_stx_rf)
print(auc_value)


# Main Figure

ggroc(list("No Space"=roc_data_stx_rf_nospace,"Spatial"=roc_data_stx_rf),
      legacy.axes = TRUE, size = 1, alpha =0.8)+
  geom_abline(linetype = "dashed") +
  scale_color_manual(values = c("Spatial" = "darkorchid2", "No Space" = "tomato1")) +  # Set custom colors
  labs(
    title = expression(paste(italic("STXBP1"), " - ROC Curve")),
    # subtitle = paste("AUC =", round(auc_value, 3)),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 14),
    text = element_text(size = 12),axis.text = element_text(size = 12)
  )

# ggsave("STXBP1_ROC_models.png",width = 6.3, height = 5, dpi = 1000)

# Test difference between spatial and non-spatial models
roc.test(roc_data_stx_rf, roc_data_stx_rf_nospace, method=c("delong"))

####################################
# SCN1A Results
####################################

scn_results = results_main_rf %>% 
  filter(gene == "SCN1A") %>% 
  mutate(actual_outcome = as.numeric(actual_outcome)) %>% 
  mutate(prediction_prob = as.numeric(prediction_prob))

scn_accuracy = sum(scn_results$hit)/nrow(scn_results)
scn_precision = sum(scn_results$TP)/(sum(scn_results$TP)+sum(scn_results$FP))
scn_recall = sum(scn_results$TP)/(sum(scn_results$TP)+sum(scn_results$FN))
scn_f1 = (2 * scn_precision * scn_recall)/(scn_precision+scn_recall)

print(paste0("Accuracy: ", scn_accuracy))

print(paste0("Precision: ",scn_precision ))

print(paste0("Recall: ",scn_recall ))

print(paste0("F1: ",scn_f1 ))

# After completing the LOOCV, compute ROC and AUC
roc_data_scn_rf <- roc(scn_results$actual_outcome, scn_results$prediction_prob)
auc_value <- auc(roc_data_scn_rf)
print(auc_value)


######
# Get ROC and AUC
######


# Main results 
# Results with localized features
scn_results = results_main_rf %>% 
  filter(gene == "SCN1A") %>% 
  mutate(actual_outcome = as.numeric(actual_outcome)) %>% 
  mutate(prediction_prob = as.numeric(prediction_prob))

roc_data_scn_rf <- roc(scn_results$actual_outcome, scn_results$prediction_prob)
auc_value <- auc(roc_data_scn_rf)
print(auc_value)

# Results with no spatial features
scn_results_nospace <- results_nospace_rf %>% 
  filter(gene == "SCN1A") %>% 
  mutate(prediction_prob = as.numeric(prediction_prob))
roc_data_scn_rf_nospace<- roc(scn_results_nospace$actual_outcome, scn_results_nospace$prediction_prob)
auc_value <- auc(roc_data_scn_rf_nospace)
print(auc_value)


ggroc(list("No Space"=roc_data_scn_rf_nospace,"Spatial"=roc_data_scn_rf),
      legacy.axes = TRUE, size = 1, alpha =0.8)+
  geom_abline(linetype = "dashed") +
  scale_color_manual(values = c("Spatial" = "#7CAE12", "No Space" = "tomato1")) +  # Set custom colors
  labs(
    title = expression(paste(italic("SCN1A"), " - ROC Curve")),
    # subtitle = paste("AUC =", round(auc_value, 3)),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 14),
    text = element_text(size = 12),axis.text = element_text(size = 12)
  )
ggsave("SCN1A_ROC_models.png",width = 6.3, height = 5, dpi = 1000)

# Test difference between ROCs
roc.test(roc_data_scn_rf, roc_data_scn_rf_nospace, method=c("delong"))

####################################
# SYNGAP1 Results
####################################

syn_results = results_main_rf %>% 
  filter(gene == "SYNGAP1") %>% 
  mutate(actual_outcome = as.numeric(actual_outcome)) %>% 
  mutate(prediction_prob = as.numeric(prediction_prob))

syn_accuracy = sum(syn_results$hit)/nrow(syn_results)
syn_precision = sum(syn_results$TP)/(sum(syn_results$TP)+sum(syn_results$FP))
syn_recall = sum(syn_results$TP)/(sum(syn_results$TP)+sum(syn_results$FN))
syn_f1 = (2 * syn_precision * syn_recall)/(syn_precision+syn_recall)

print(paste0("Accuracy: ", syn_accuracy))

print(paste0("Precision: ",syn_precision ))

print(paste0("Recall: ",syn_recall ))

print(paste0("F1: ",syn_f1 ))

######
# ROC and AUC Results
######


# Main results
# Results with spatial features
syn_results = results_main_rf %>% 
  filter(gene == "SYNGAP1") %>% 
  mutate(actual_outcome = as.numeric(actual_outcome)) %>% 
  mutate(prediction_prob = as.numeric(prediction_prob))

roc_data_syn_rf <- roc(syn_results$actual_outcome, syn_results$prediction_prob)
auc_value <- auc(roc_data_syn_rf)
print(auc_value)

# Results with no spatial features
syn_results_nospace <- results_nospace_rf %>% 
  filter(gene == "SYNGAP1") %>% 
  mutate(prediction_prob = as.numeric(prediction_prob))
roc_data_syn_rf_nospace<- roc(syn_results_nospace$actual_outcome, syn_results_nospace$prediction_prob)
auc_value <- auc(roc_data_syn_rf_nospace)
print(auc_value)


ggroc(list("No Space"=roc_data_syn_rf_nospace,"Spatial"=roc_data_syn_rf),
      legacy.axes = TRUE, size = 1, alpha =0.8)+
  geom_abline(linetype = "dashed") +
  scale_color_manual(values = c("Spatial" = "#00BFC9", "No Space" = "tomato1")) +  # Set custom colors
  labs(
    title = expression(paste(italic("SYNGAP1"), " - ROC Curve")),
    # subtitle = paste("AUC =", round(auc_value, 3)),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 14),
    text = element_text(size = 12),axis.text = element_text(size = 12)
  )
ggsave("SYNGAP1_ROC_models.png",width = 6.3, height = 5, dpi = 1000)

# Test difference between ROCs
roc.test(roc_data_syn_rf, roc_data_syn_rf_nospace, method=c("delong"))


########################################################
# Top Features
########################################################

###############
# For just lobe features models
###############

results_long_pre <- results_main_rf %>%
  mutate(across(.cols = -gene, .fns = ~ as.numeric(as.character(.)))) %>% 
  mutate(age = as.numeric(age)) %>%
  select(-c(hit,FP, TP, TN, FN)) %>% 
  pivot_longer(
    cols = -c(1:6), # Keep the first 4 columns fixed; adjust the number if needed
    names_to = "features", # Name of the column to store the original column names
    values_to = "weights"  # Name of the column to store the values
  ) %>% 
  mutate(features = gsub("Beta Theta","Beta-Theta",features)) %>% 
  mutate(features = gsub("Alpha Delta","Alpha-Delta",features)) %>% 
  mutate(features = gsub("Alpha Theta","Alpha-Theta",features)) 


results_long_max <- results_long_pre %>% 
  group_by(gene, Iteration) %>%
  # To aid in normalization
  dplyr::summarise(
    min_weight = min(weights),
    max_weight = max(weights)) %>% 
  ungroup()


results_long<- results_main_rf %>%
  mutate(across(.cols = -gene, .fns = ~ as.numeric(as.character(.)))) %>% 
  mutate(age = as.numeric(age)) %>%
  select(-c(hit,FP, TP, TN, FN)) %>% 
  pivot_longer(
    cols = -c(1:5), # Keep the first 4 columns fixed; adjust the number if needed
    names_to = "features", # Name of the column to store the original column names
    values_to = "weights"  # Name of the column to store the values
  ) %>% 
  left_join(results_long_max) %>% 
  # For mean reduced accuracy
  mutate(shifted_importance = weights - abs(min_weight)) %>%  # Shift all values
  group_by(gene, Iteration) %>%
  dplyr::summarise(
    max_shifted_weight = max(shifted_importance),
    min_shifted_weight = min(shifted_importance)) %>% 
  right_join(results_long_pre) %>% 
  left_join(results_long_max) %>% 
  mutate(shifted_importance = weights - abs(min_weight)) %>%  # Shift all values
  mutate(normalized_weight = (shifted_importance - min_shifted_weight) / (max_shifted_weight - min_shifted_weight))



###############
# For just all features models
###############
results_long_pre <- results_all_rf %>%
  mutate(across(.cols = -gene, .fns = ~ as.numeric(as.character(.)))) %>% 
  mutate(age = as.numeric(age)) %>%
  select(-c(hit,FP, TP, TN, FN)) %>% 
  pivot_longer(
    cols = -c(1:6), # Keep the first 4 columns fixed; adjust the number if needed
    names_to = "features", # Name of the column to store the original column names
    values_to = "weights"  # Name of the column to store the values
  ) %>% 
  mutate(features = gsub("Beta Theta","Beta-Theta",features)) %>% 
  mutate(features = gsub("Alpha Delta","Alpha-Delta",features)) %>% 
  mutate(features = gsub("Alpha Theta","Alpha-Theta",features)) 


results_long_max <- results_long_pre %>% 
  group_by(gene, Iteration) %>%
  # For mean reduced accuracy
  dplyr::summarise(
    min_weight = min(weights),
    max_weight = max(weights)) %>% 
  ungroup()


results_all_long<- results_all_rf %>%
  mutate(across(.cols = -gene, .fns = ~ as.numeric(as.character(.)))) %>% 
  mutate(age = as.numeric(age)) %>%
  select(-c(hit,FP, TP, TN, FN)) %>% 
  pivot_longer(
    cols = -c(1:5), # Keep the first 4 columns fixed; adjust the number if needed
    names_to = "features", # Name of the column to store the original column names
    values_to = "weights"  # Name of the column to store the values
  ) %>% 
  left_join(results_long_max) %>% 
  # For mean reduced accuracy
  mutate(shifted_importance = weights - abs(min_weight)) %>%  # Shift all values
  group_by(gene, Iteration) %>%
  dplyr::summarise(
    max_shifted_weight = max(shifted_importance),
    min_shifted_weight = min(shifted_importance)) %>% 
  right_join(results_long_pre) %>% 
  left_join(results_long_max) %>% 
  mutate(shifted_importance = weights - abs(min_weight)) %>%  # Shift all values
  mutate(normalized_weight = (shifted_importance - min_shifted_weight) / (max_shifted_weight - min_shifted_weight))



############
# STXBP1
############

# Just Lobe Features
top_stx <- results_long %>%
  filter(gene == "STXBP1") %>% 
  group_by(gene, features) %>% 
  dplyr::summarise(med_weight = median(normalized_weight)) %>% 
  arrange(desc(abs(med_weight))) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  select(gene, med_weight,features) %>%
  mutate(top = "Yes")

stx_features  <- results_long %>%
  filter(gene == "STXBP1") %>% 
  left_join(top_stx) %>% 
  filter(top == "Yes") %>% 
  mutate(features = gsub("_"," ", features)) %>% 
  mutate(features = str_to_title(features)) %>% 
  mutate(features = gsub("Beta Theta","Beta-Theta",features)) %>% 
  mutate(features = gsub("Alpha Delta","Alpha-Delta",features)) %>% 
  mutate(features = gsub("Beta Delta","Beta-Delta",features)) %>% 
  
  mutate(features = gsub("Alpha Theta","Alpha-Theta",features)) 


# Sort by median feature weight
ordered_features <- stx_features %>%
  select(features, med_weight) %>% 
  distinct() %>% 
  group_by(features) %>% 
  arrange(desc(med_weight)) %>%
  pull(features)

# Reorder the feature factor in the entire dataset based on the above ordering
stx_features$features <- factor(stx_features$features, levels = ordered_features)


ggplot(stx_features, aes(features, normalized_weight))+
  geom_boxplot(color = "#C77CFF", fill = "#C77CFF",alpha=0.25)+
  ylab("Normalized Feature Weight")+
  ggtitle(expression(italic("STXBP1") ~ "Random Forest Feature Weights"))+
  scale_y_continuous(breaks = seq(0.2, 1, len = 5)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_blank()  # Remove x-axis title
  )

ggsave("STXBP1_RF_model_weights_GINI.png",width = 6, height = 5, dpi = 1000)


################
# SCN1A
################


# Just Lobe Features

top_scn <- results_long %>%
  filter(gene == "SCN1A") %>% 
  group_by(gene, features) %>% 
  dplyr::summarise(med_weight = median(normalized_weight)) %>% 
  arrange(desc(abs(med_weight))) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  select(gene, med_weight,features) %>%
  mutate(top = "Yes")

scn_features  <- results_long %>%
  filter(gene == "SCN1A") %>% 
  left_join(top_scn) %>% 
  filter(top == "Yes") %>% 
  mutate(features = gsub("_"," ", features)) %>% 
  mutate(features = str_to_title(features))%>% 
  mutate(features = gsub("Beta Theta","Beta-Theta",features)) %>% 
  mutate(features = gsub("Alpha Delta","Alpha-Delta",features)) %>% 
  mutate(features = gsub("Beta Delta","Beta-Delta",features)) %>% 
  
  mutate(features = gsub("Alpha Theta","Alpha-Theta",features)) 


# Filter to "Infant" age group and arrange by kw_chi
ordered_features <- scn_features %>%
  select(features, med_weight) %>% 
  distinct() %>% 
  group_by(features) %>% 
  arrange(desc(med_weight)) %>%
  pull(features)

# Reorder the feature factor in the entire dataset based on the above ordering
scn_features$features <- factor(scn_features$features, levels = ordered_features)


ggplot(scn_features, aes(features, normalized_weight))+
  geom_boxplot(color = "#7CAE00", fill = "#7CAE00",alpha=0.25)+
  ylab("Normalized Feature Weight")+
  ggtitle(expression(italic("SCN1A") ~ "Random Forest Feature Weights"))+
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_blank()  # Remove x-axis title
  )

ggsave("SCN1A_RF_model_weights_GINI.png",width = 6, height = 5, dpi = 1000)

############
# SYNGAP1
############

# Just Lobe Features
top_syn <- results_long %>%
  filter(gene == "SYNGAP1") %>% 
  group_by(gene, features) %>% 
  dplyr::summarise(med_weight = median(normalized_weight)) %>% 
  arrange(desc(abs(med_weight))) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  select(gene, med_weight,features) %>%
  mutate(top = "Yes")

syn_features  <- results_long %>%
  filter(gene == "SYNGAP1") %>% 
  left_join(top_syn) %>% 
  filter(top == "Yes") %>% 
  mutate(features = gsub("_"," ", features)) %>% 
  mutate(features = str_to_title(features)) %>% 
  mutate(features = gsub("Beta Theta","Beta-Theta",features)) %>% 
  mutate(features = gsub("Alpha Delta","Alpha-Delta",features)) %>% 
  mutate(features = gsub("Beta Delta","Beta-Delta",features)) %>% 
  
  mutate(features = gsub("Alpha Theta","Alpha-Theta",features)) 

# Filter to "Infant" age group and arrange by kw_chi
ordered_features <- syn_features %>%
  select(features, med_weight) %>% 
  distinct() %>% 
  group_by(features) %>% 
  arrange(desc(med_weight)) %>%
  pull(features)

# Reorder the feature factor in the entire dataset based on the above ordering
syn_features$features <- factor(syn_features$features, levels = ordered_features)


ggplot(syn_features, aes(features, normalized_weight))+
  geom_boxplot(color = "#00BFC4", fill = "#00BFC4", alpha=0.25)+
  ylab("Normalized Feature Weight")+
  ggtitle(expression(italic("SYNGAP1") ~ "Random Forest Feature Weights"))+
  scale_y_continuous(breaks = round(seq(0.2, max(syn_features$normalized_weight)+0.53, by = 0.2),1)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_blank()  # Remove x-axis title
  )

ggsave("SYNGAP1_RF_model_weights_GINI.png",width = 6, height = 5, dpi = 1000)


#################
# 3-WAY comparison with resampling and bootstrapping
#################  

# Assuming feat_lobe_wide is your dataframe and has been prepared up to this point:
data_pretrain <- feat_lobe_wide %>%
  filter(gene %in% c("SCN1A", "SYNGAP1", "STXBP1")) %>% 
  mutate(gene = as.factor(gene))

# Extract data for each gene
scn1a_data <- data_pretrain %>% filter(gene == "SCN1A")
syngap1_data <- data_pretrain %>% filter(gene == "SYNGAP1")
stxbp1_data <- data_pretrain %>% filter(gene == "STXBP1")

# Number of samples equal to the number of SYNGAP1 samples
n_samples <- nrow(syngap1_data)

# Initialize dataframe to store results
results_df <- as.data.frame(matrix(nrow=0, ncol=6))
names(results_df) <- c("Iteration", "prediction_prob", "prediction", 
                       "actual_outcome", "age_test", "age")

# Sample with replacement and apply 80/20 split for each gene
sample_and_split <- function(data) {
  sample_indices <- sample(1:nrow(data), size = n_samples, replace = TRUE)
  split_index <- floor(0.8 * n_samples)
  list(
    train = data[sample_indices[1:split_index], ],
    test = data[sample_indices[(split_index + 1):n_samples], ]
  )
}
set.seed(123)  # For reproducibility

# To store results from each iteration
results_df <- list()
importance_df <- list()

#Set number of bootstraps
n_bootstraps <- 10000

for (i in 1:n_bootstraps) {
  
  scn1a_split <- sample_and_split(scn1a_data)
  syngap1_split <- sample_and_split(syngap1_data)
  stxbp1_split <- sample_and_split(stxbp1_data)
  
  # Combine training data from each gene
  training_data <- bind_rows(scn1a_split$train, syngap1_split$train, stxbp1_split$train)
  testing_data <- bind_rows(scn1a_split$test, syngap1_split$test, stxbp1_split$test)
  
  # Train Random Forest model
  rf_model <- randomForest(gene ~ ., data = training_data %>% select(-c(patient, age_group)), ntree = 500)
  importance_vals <- importance(rf_model) %>% as.data.frame() %>% t()
  rownames(importance_vals) <- NULL
  
  # Train Null Age Model
  null_model <- randomForest(gene ~ ., data = training_data %>% select(gene, age), ntree = 500)
  
  ############
  # Predict and evaluate
  predict_prob <- predict(rf_model, newdata = testing_data,type = "prob")
  prediction <- predict(rf_model, newdata = testing_data)
  # Calculate ROC and AUC for each class against all others
  roc_results <- lapply(levels(testing_data$gene), function(class) {
    binary_labels <- ifelse(testing_data$gene == class, class, "Others")
    roc_curve <- roc(response = binary_labels, predictor = predict_prob[, class])
    auc_value <- auc(roc_curve)
    list(roc_curve = roc_curve, auc = auc_value)
  })
  
  stxbp1_auc <- roc_results[[which(levels(testing_data$gene) == "STXBP1")]]$auc[1]
  syngap1_auc <- roc_results[[which(levels(testing_data$gene) == "SYNGAP1")]]$auc[1]
  scn1a_auc <- roc_results[[which(levels(testing_data$gene) == "SCN1A")]]$auc[1]
  ############
  ############
  null_predict_prob <- predict(null_model, newdata = testing_data,type = "prob")
  null_prediction <- predict(null_model, newdata = testing_data)
  # Calculate ROC and AUC for each class against all others
  roc_results <- lapply(levels(testing_data$gene), function(class) {
    binary_labels <- ifelse(testing_data$gene == class, class, "Others")
    roc_curve <- roc(response = binary_labels, predictor = null_predict_prob[, class])
    auc_value <- auc(roc_curve)
    list(roc_curve = roc_curve, auc = auc_value)
  })
  
  # Get Null AUCs
  null_stxbp1_auc <- roc_results[[which(levels(testing_data$gene) == "STXBP1")]]$auc[1]
  null_syngap1_auc <- roc_results[[which(levels(testing_data$gene) == "SYNGAP1")]]$auc[1]
  null_scn1a_auc <- roc_results[[which(levels(testing_data$gene) == "SCN1A")]]$auc[1]
  #######################
  
  actual <- testing_data$gene
  
  # Collect results
  cur_df <- data.frame(
    Iteration = i,
    prediction = prediction,
    null_prediction = null_prediction,
    actual_outcome = actual,
    stx_auc = stxbp1_auc,
    syn_auc = syngap1_auc,
    scn_auc = scn1a_auc,
    null_stx_auc = null_stxbp1_auc,
    null_syn_auc = null_syngap1_auc,
    null_scn_auc = null_scn1a_auc,
    age_test = testing_data$age  
  )
  results_df[[i]] <- cur_df
  importance_df[[i]] <- importance_vals %>% as.data.frame()
}

# Combine all results into a single data frame
final_results_df <- bind_rows(results_df)
final_importance_df <- bind_rows(importance_df)

###############
# Overall Accuracy
###############
# Calculate accuracy for each iteration
accuracy_results <- final_results_df %>%
  group_by(Iteration) %>%
  summarise(
    main_model_acc = mean(prediction == actual_outcome),  # Accuracy of the primary model
    null_model_acc = mean(null_prediction == actual_outcome)  # Accuracy of the null model
  )

# Calculate the median accuracy across all bootstraps for primary model and null model
median_main_accuracy <- median(accuracy_results$main_model_acc)
median_null_accuracy <- median(accuracy_results$null_model_acc)

# Print the median accuracy 
print(paste("Median accuracy for primary model:", median_main_accuracy))
print(paste("Median accuracy for null model:", median_null_accuracy))


###############
# AUC Values
###############
summary(final_results_df$null_stx_auc)
summary(final_results_df$null_scn_auc)
summary(final_results_df$null_syn_auc)


summary(final_results_df$stx_auc)
summary(final_results_df$scn_auc)
summary(final_results_df$syn_auc)


basic_conf_matrix <- table(Predicted = final_results_df$prediction, Actual = final_results_df$actual_outcome)
print(basic_conf_matrix)

# Ensure that both predicted and actual outcome columns are factors and have the same levels
final_results_df$prediction <- factor(final_results_df$prediction, levels = unique(final_results_df$actual_outcome))
final_results_df$actual_outcome <- factor(final_results_df$actual_outcome)

# Generate confusion matrix
detailed_conf_matrix <- confusionMatrix(data = final_results_df$prediction, reference = final_results_df$actual_outcome)
print(detailed_conf_matrix)

# Extract the byClass statistics
class_stats <- detailed_conf_matrix$byClass

# Calculate F1 score for each class
f1_scores <- 2 * (class_stats[,"Sensitivity"] * class_stats[,"Precision"]) / (class_stats[,"Sensitivity"] + class_stats[,"Precision"])


###############
# Confusion Matrix Plot
###############

results_list <- final_results_df %>%
  mutate(
    prediction = factor(prediction, levels = c("STXBP1", "SYNGAP1", "SCN1A")),
    actual_outcome = factor(actual_outcome, levels = c("SCN1A", "SYNGAP1", "STXBP1"))
  ) %>% 
  group_by(Iteration) %>%
  summarise(conf_mat = list(prop.table(table(prediction, actual_outcome), margin = 2))) %>%
  pull(conf_mat)

# Summarize all percentage confusion matrices to calculate mean percentage matrix
total_conf_mat <- Reduce(`+`, results_list)
mean_percentage_conf_mat <- total_conf_mat / length(results_list)

# Calculate the standard deviation for each cell across all percentage confusion matrices
sd_mat <- Reduce("+", lapply(results_list, function(x) (x - mean_percentage_conf_mat)^2)) / (length(results_list) - 1)

se_mat <- sqrt(sd_mat) / sqrt(length(results_list))

# 95% Confidence Interval 
ci_lower <- mean_percentage_conf_mat - 1.96 * se_mat
ci_upper <- mean_percentage_conf_mat + 1.96 * se_mat

# Assuming mean_percentage_conf_mat, ci_lower, and ci_upper are matrices
mean_conf_df <- as.data.frame(mean_percentage_conf_mat)
names(mean_conf_df) <- c("Prediction", "Actual", "Mean")

ci_lower_df <- as.data.frame(ci_lower, row.names = NULL)
ci_upper_df <- as.data.frame(ci_upper, row.names = NULL)
names(ci_lower_df) <- names(ci_upper_df) <- c("Prediction", "Actual", "CI")

# Merging confidence intervals with the mean data frame
mean_conf_df$CI_lower <- ci_lower_df$CI
mean_conf_df$CI_upper <- ci_upper_df$CI

# Flatten the data frame for ggplot (if matrices were reshaped incorrectly)
mean_conf_df <- reshape2::melt(mean_percentage_conf_mat)
names(mean_conf_df) <- c("Actual", "Prediction", "Mean")
mean_conf_df$CI_lower <- melt(ci_lower)$value
mean_conf_df$CI_upper <- melt(ci_upper)$value

mean_conf_df <- mean_conf_df %>% 
  mutate(
    Prediction = factor(Prediction, levels = c("STXBP1", "SYNGAP1", "SCN1A")),
    Actual = factor(Actual, levels = c("SCN1A", "SYNGAP1", "STXBP1")))

# Plotting the confusion matrix
ggplot(mean_conf_df, aes(Prediction, Actual, fill = Mean)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f%%\n[%.2f%%, %.2f%%]", Mean*100, CI_lower*100, CI_upper*100)), size = 3) +
  scale_fill_gradient(low = "white", high = "dodgerblue2") +
  labs(title = "Average Confusion Matrix with 95% CI", x = "Predicted Class", y = "Actual Class") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_blank()  # Optionally remove x-axis title
  )

# ggsave("All_gene_3x_RF_model_conf_matrix_labels_perc.png",width = 6.3, height = 5, dpi = 1000)

####################
# Feature Importance
####################

final_importance_df$Iteration <- seq(1, 1000, 1)
importance_max <- final_importance_df %>% 
  rowwise() %>%  # Apply functions to each row
  mutate(
    min_weight = min(c_across(-Iteration), na.rm = TRUE),  # Exclude 'Iteration' column
    max_weight = max(c_across(-Iteration), na.rm = TRUE)   # Exclude 'Iteration' column
  ) %>%
  ungroup() %>%  # Remove rowwise grouping for further data manipulation
  select(Iteration, min_weight, max_weight)

importance_long<- final_importance_df %>%
  
  pivot_longer(
    cols = -c(Iteration), # Keep the first 4 columns fixed; adjust the number if needed
    names_to = "features", # Name of the column to store the original column names
    values_to = "weights"  # Name of the column to store the values
  ) 

importance_long_norm <- importance_long %>% 
  left_join(importance_max) %>% 
  # For mean reduced accuracy
  mutate(shifted_importance = weights - abs(min_weight)) %>%  # Shift all values
  group_by(Iteration) %>%
  dplyr::summarise(
    max_shifted_weight = max(shifted_importance),
    min_shifted_weight = min(shifted_importance)) %>% 
  right_join(importance_long) %>%
  left_join(importance_max) %>% 
  mutate(shifted_importance = weights - abs(min_weight)) %>%  # Shift all values
  mutate(normalized_weight = (shifted_importance - min_shifted_weight) / (max_shifted_weight - min_shifted_weight))


# ALL Features
top_feats <- importance_long_norm %>%
  group_by(features) %>% 
  dplyr::summarise(med_weight = median(normalized_weight)) %>% 
  arrange(desc(abs(med_weight))) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  select(med_weight,features) %>%
  mutate(top = "Yes")


all_features  <- importance_long_norm %>%
  left_join(top_feats) %>% 
  filter(top == "Yes") %>% 
  mutate(features = gsub("_"," ", features)) %>% 
  mutate(features = str_to_title(features)) %>% 
  mutate(features = gsub("Beta Theta","Beta-Theta",features)) %>% 
  mutate(features = gsub("Alpha Delta","Alpha-Delta",features)) %>% 
  mutate(features = gsub("Beta Delta","Beta-Delta",features)) %>% 
  
  mutate(features = gsub("Alpha Theta","Alpha-Theta",features)) 

# Filter to "Infant" age group and arrange by kw_chi
ordered_features <- all_features %>%
  select(features, med_weight) %>% 
  distinct() %>% 
  group_by(features) %>% 
  arrange(desc(med_weight)) %>%
  pull(features)

# Reorder the feature factor in the entire dataset based on the above ordering
all_features$features <- factor(all_features$features, levels = ordered_features)


ggplot(all_features, aes(features, normalized_weight))+
  geom_boxplot(color = "dodgerblue2", fill = "dodgerblue2", alpha=0.25)+
  ylab("Normalized Feature Weight")+
  ggtitle("3-Way Random Forest Feature Weights")+
  scale_y_continuous(breaks = round(seq(0.2, max(all_features$normalized_weight)+0.53, by = 0.2),1)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_blank()  # Remove x-axis title
  )

# ggsave("All_gene_3x_RF_model_weights_all_feats_gini.png",width = 6, height = 5, dpi = 1000)

