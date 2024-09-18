####
# Description: 
# This file creates the random forest models to predict seizure frequency. 
# These results are seen in Figure 5 of Galer et al.
# Note, results may differ very slightly due to the slight randomness of random 
# forest models
####

# Loads primary functions
source("scripts/main_R_functions.R")



#### Load Data ####
# Spectral information for every electrode

# Create wide format, remove unnecessary columns, convert outcome to binary
# Note all of the patients have a known or presumed genetic epilepsy roughly defined here as "DEE"
full_seiz_psd <- read_csv("data/psd_bp_seiz_freq.csv") %>% 
  pivot_wider(names_from = feature, values_from = med_feat,
              names_prefix = "", id_cols = c(patient, age, seiz_age,seiz_freqs)) %>% 
  # This can be manipulated based on how you want to split categories
  mutate(broad_freq = case_when((seiz_freqs %in% c(">1 per day",">1 per week")) ~ "Frequent",
                                (seiz_freqs %in% c( ">2 years ago","1-2 years ago",">1 per 6 months",'>1 per 1 year',">1 per month"))~"Rare",
                                TRUE~"Other")) %>% 
  filter(broad_freq !="Other") %>% 
  mutate(seiz_binary = ifelse(broad_freq == "Frequent", 1, 0)) %>% 
  mutate(seiz_binary = as.factor(seiz_binary))



##### Histogram -----

ordered_features <- c(">1 per day",">1 per week",">1 per month",">1 per 6 months",'>1 per 1 year',"1-2 years ago", ">2 years ago")


# Look at different frequencies
full_seiz_psd$seiz_freqs <- factor(full_seiz_psd$seiz_freqs, levels = ordered_features)



#### All DEEs ####
seiz_wide_cnt = full_seiz_psd %>% 
  count(seiz_freqs)
seiz_wide_cnt$seiz_freqs <- factor(seiz_wide_cnt$seiz_freqs, levels = ordered_features)

ggplot(seiz_wide_cnt, aes(seiz_freqs, n))+
  geom_histogram(stat = "identity", fill = "maroon", alpha=0.9) +
  ggtitle("DEE Seizure Frequency Scores")+
  scale_y_continuous(breaks = round(seq(0, max(seiz_wide_cnt$n), by = 25),0)) +
  ylab("Count") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_blank()  # Remove x-axis title
  )

ggsave("seiz_hist_DEE_pats.png",width = 5, height = 5, dpi = 1000)

#### Random Forest Model ####


train_sz_freq <- full_seiz_psd %>% 
  select(-broad_freq)

# Determine which is larger - larger becomes "control"
## This is how the find_age_matched_control() works
length(unique(train_sz_freq$patient[train_sz_freq$seiz_binary==0]))
length(unique(train_sz_freq$patient[train_sz_freq$seiz_binary==1]))

# This keeps the same syntax as other analyses (gene vs control)
df_gene <- train_sz_freq %>% 
  filter(seiz_binary==1)
df_control <- train_sz_freq %>% 
  filter(seiz_binary==0)


# Perform age matching
control_match <- find_age_matched_controls(df_gene, df_control)

# Check to make sure that no patient is duplicated
if (length(unique(control_match$patient)) != nrow(control_match)) {
  print("WARNING: One or more patients are duplicated. Check code.")
}


# Combine "Controls" with "Gene"
train_sz_freq <- df_gene %>% 
  rbind(control_match)

#### Performing LOOCV ####
sz_freq_results <- seizfreq_rf_loocv(train_sz_freq)


sz_freq_results_rf = sz_freq_results %>% 
  mutate(Iteration = as.numeric(Iteration)) %>% 
  mutate(actual_outcome = (as.numeric(actual_outcome)-1)) %>%
  mutate(prediction = as.numeric(prediction)) %>%
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

sz_freq_results_rf_null = sz_freq_results %>% 
  mutate(Iteration = as.numeric(Iteration)) %>% 
  mutate(actual_outcome = (as.numeric(actual_outcome)-1)) %>%
  mutate(prediction = as.numeric(null_prediction)) %>%
  mutate(prediction_prob = as.numeric(null_prediction_prob)) %>% 
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

sz_freq_accuracy = sum(sz_freq_results_rf$hit)/nrow(sz_freq_results_rf)
sz_freq_precision = sum(sz_freq_results_rf$TP)/(sum(sz_freq_results_rf$TP)+sum(sz_freq_results_rf$FP))
sz_freq_recall = sum(sz_freq_results_rf$TP)/(sum(sz_freq_results_rf$TP)+sum(sz_freq_results_rf$FN))
sz_freq_f1 = (2 * sz_freq_precision * sz_freq_recall)/(sz_freq_precision+sz_freq_recall)

sz_freq_accuracy_null = sum(sz_freq_results_rf_null$hit)/nrow(sz_freq_results_rf_null)
sz_freq_precision_null = sum(sz_freq_results_rf_null$TP)/(sum(sz_freq_results_rf_null$TP)+sum(sz_freq_results_rf_null$FP))
sz_freq_recall_null = sum(sz_freq_results_rf_null$TP)/(sum(sz_freq_results_rf_null$TP)+sum(sz_freq_results_rf_null$FN))
sz_freq_f1_null = (2 * sz_freq_precision_null * sz_freq_recall_null)/(sz_freq_precision_null+sz_freq_recall_null)

print(paste0("Accuracy: ", sz_freq_accuracy))
print(paste0("Null Accuracy: ", sz_freq_accuracy_null))

print(paste0("Precision: ",sz_freq_precision ))
print(paste0("Null Precision: ",sz_freq_precision_null))

print(paste0("Recall: ",sz_freq_recall ))
print(paste0("Null Recall: ",sz_freq_recall_null ))

print(paste0("F1: ",sz_freq_f1 ))
print(paste0("Null F1: ",sz_freq_f1_null ))

roc_data_rf <- roc(sz_freq_results_rf$actual_outcome, sz_freq_results_rf$prediction_prob)
auc_value <- auc(roc_data_rf)
print(auc_value)


roc_data_rf_null <- roc(sz_freq_results_rf_null$actual_outcome, sz_freq_results_rf_null$null_prediction_prob)
auc_value <- auc(roc_data_rf_null)
print(auc_value)

ggroc(roc_data_rf,
      legacy.axes = TRUE, size = 1, alpha =1, color = "maroon")+
  geom_abline(linetype = "dashed") +
  # scale_color_manual(values = c("Spatial" = "darkorchid2", "No Space" = "tomato1")) +  # Set custom colors
  labs(
    title = ("Seizure Frequency - ROC Curve"),
    # subtitle = paste("AUC =", round(auc_value, 3)),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16),
    plot.subtitle = element_text(size = 14),
    text = element_text(size = 12),axis.text = element_text(size = 12)
  )

ggsave("seiz_rf_loocv_ROC_sz_freq.png",width = 5.1, height = 5, dpi = 1000)
