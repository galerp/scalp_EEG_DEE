###################
# Description: 
# This file creates the random forest models to predict Gross Motor Measure 
# (GMFM) scores. These results are seen in Figure 5 of Galer et al.
# Note, results may differ very slightly due to randomness of random forest models
# and random sampling each bootstrap. e.g., the top 8 features may shift 
# slightly in order (ordered by median), but should stay the same.

###################

# Set working directory
setwd("/Volumes/helbig_lab/Users/galerp/EEG/manuscript/github/")

# Loads primary functions
source("main_R_functions.R")

# Spectral information for every electrode

psd_bp_rel_gmfm <-  read_csv("data/psd_bp_GMFM.csv")

# Find median for each feature at each lobe
psd_bp_rel_lobe <- psd_bp_rel_gmfm %>% 
  filter(electrode %nin% c("A1","A2")) %>% 
  group_by(patient, MEAS_VALUE, age_test, age_eeg, age_group,lobe, bandpower, epochs) %>% 
  dplyr::summarise(med_feat = median(power)) %>% 
  mutate(feature = paste(bandpower, lobe, sep = "_")) %>% 
  mutate(feature = gsub("_power_rel","",feature)) %>% 
  mutate(feature = gsub("Frontal","Frnt",feature)) %>% 
  mutate(feature = gsub("Occipital","Occ",feature)) %>% 
  mutate(feature = gsub("Temporal","Tmp",feature)) %>% 
  mutate(feature = gsub("Central","Cnt",feature)) %>% 
  mutate(feature = gsub("Parietal","Par",feature)) %>% 
  select(-lobe) 


# Convert to wide format
psd_rel_rat_wide <- psd_bp_rel_lobe %>%
  pivot_wider(names_from = feature, values_from = med_feat,
              names_prefix = "", id_cols = c(patient, MEAS_VALUE,age_test,age_eeg, age_group, epochs))

######################
# Train and Test
######################
train_df <- psd_rel_rat_wide %>%
  filter(epochs>=15) %>% 
  select(-c(patient, epochs, age_eeg, age_group)) %>% 
  ungroup()

# Note this can take 10-20 minutes to run
rf_results <- gmfm_rf(train_df)

# Take important results
rf_comp <- rf_results %>% 
  select(iteration, rmse, null_rmse, rmse_dif, mae, null_mae, mae_dif,dif_comp) %>% 
  distinct()

# Find if EEG models performs significantly better via:
# MAE - mean absolued error
t.test(rf_comp$mae_dif,rf_comp$dif_comp,alternative = "greater")
# RMSE - root mean squared error (primrary measure in paper)
t.test(rf_comp$rmse_dif,rf_comp$dif_comp,alternative = "greater")

# Calculate additional statistics
# median differences
median_mae_dif <- median(rf_comp$mae_dif)
median_rmse_diff <- median(rf_comp$rmse_dif)
# median null performance
null_mae <- median(rf_comp$null_mae)
null_rmse <- median(rf_comp$null_rmse)
# median EEG performance
eeg_mae <- median(rf_comp$mae)
eeg_rmse <- median(rf_comp$rmse)


#############
# Plot RMSE Difference
# As seen in Figure 5
#############


rf_rmse_long <- rf_comp %>% 
  select(iteration, rmse) %>% 
  mutate(Model = "EEG") %>% 
  rbind(rf_comp %>% 
          select(iteration, null_rmse) %>% 
          rename(rmse = null_rmse) %>% 
          mutate(Model = "Null"))
rf_rmse_long$Model = factor(rf_rmse_long$Model, levels = c("Null","EEG"))


# Find the median difference
median_difference <- median(rf_comp$rmse_dif)

# Find the iteration with the difference closest to the median difference
closest_iteration <- rf_comp %>%
  mutate(Closest_to_Median = abs(rmse_dif - median_difference)) %>%
  dplyr::arrange(Closest_to_Median) %>%
  slice(1) %>%
  pull(iteration)


ggplot(rf_rmse_long, aes(x = Model, y = rmse, group = iteration)) +
  geom_point(color = "cadetblue3", position = position_dodge(width = 0.25), size = 1, alpha=0.3) +
  geom_line( color = "cadetblue3",position = position_dodge(width = 0.25), alpha = 0.03) +
  # theme_minimal() + 
  geom_line(data = rf_rmse_long %>% filter(iteration == closest_iteration),
            aes(x = Model, y = rmse),
            color = "red3", size = 1.5, alpha = 0.8)+
  xlab("Model")+
  ylab("RMSE")+
  # ylim(c(0,100))+
  scale_x_discrete(expand = c(0.15, 0.15)) +
  ylim(c(11.5,33))+
  ggtitle("GMFM Model RMSE") +
  theme_classic()+
  theme(text = element_text(size = 14),axis.text = element_text(size = 14))


ggsave("All_gmfm_model_RMSE_performance.png",
       width = 3.2, height = 5, dpi=1000)


#########
# Feature Importance
#########

results_long_pre <- rf_results %>%
  select(-c(prediction, null_predictions, actual_outcome, age, mae, null_mae, mae_dif, rmse, null_rmse, rmse_dif, dif_comp, age_test)) %>%
  pivot_longer(
    cols = -c(1), 
    names_to = "features", # Name of the column to store the original column names
    values_to = "weights"  # Name of the column to store the values
  ) %>% 
  distinct()


results_long_max <- results_long_pre %>% 
  group_by(iteration) %>%
  # For mean reduced accuracy
  dplyr::summarise(
    min_weight = min(weights),
    max_weight = max(weights)) %>% 
  ungroup()
range(results_long_max$min_weight)

results_long<- results_long_pre %>%
  left_join(results_long_max) %>% 
  # For mean reduced accuracy
  mutate(shifted_importance = weights - min_weight) %>%  
  group_by(iteration) %>%
  dplyr::summarise(
    max_shifted_weight = max(shifted_importance),
    min_shifted_weight = min(shifted_importance)) %>% 
  right_join(results_long_pre) %>% 
  left_join(results_long_max) %>% 
  mutate(shifted_importance = weights - min_weight) %>%  
  mutate(normalized_weight = (shifted_importance - min_shifted_weight) / (max_shifted_weight - min_shifted_weight)) %>% 
  mutate(features = gsub("_"," ",features)) %>% 
  mutate(features = str_to_title(features)) %>% 
  mutate(features = gsub("Beta Theta","Beta-Theta",features)) %>% 
  mutate(features = gsub("Alpha Delta","Alpha-Delta",features)) %>% 
  mutate(features = gsub("Alpha Theta","Alpha-Theta",features)) 

top_feats <- results_long %>%
  group_by(features) %>% 
  dplyr::summarise(med_weight = median(normalized_weight)) %>% 
  dplyr::arrange(desc(abs(med_weight))) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  select(med_weight,features) %>%
  mutate(top = "Yes")


# To get correct order (sorted by median)
ordered_features <- top_feats %>%
  select(features, med_weight) %>% 
  distinct() %>% 
  group_by(features) %>% 
  dplyr::arrange(desc(med_weight)) %>%
  pull(features)

top_feats_weights <- results_long %>% 
  filter(features%in% ordered_features)

# Reorder the feature factor based on the above ordering
top_feats_weights$features <- factor(top_feats_weights$features, levels = ordered_features)


ggplot(top_feats_weights, aes(features, normalized_weight))+
  geom_boxplot(color = "cadetblue3", fill = "cadetblue3", alpha = 0.5)+
  ylab("Normalized Feature Weight")+
  ggtitle("GMFM Random Forest Feature Weights")+
  scale_y_continuous(limits = c(0,1))+
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_blank()  # Remove x-axis title
  )



# ggsave("All_gmfm_model_feat_importance_GINI.png",width = 5.9, height = 5.3, dpi = 1000)




####################
# More statistics if wanted
####################
# Compute size
n <- nrow(rf_comp)
# Find the standard deviation
mae_sd <- sd(rf_comp$mae_dif)
rmse_sd <- sd(rf_comp$rmse_dif)

# Find the standard error
se_mae <- mae_sd / sqrt(n)
se_rmse <- rmse_sd / sqrt(n)

alpha = 0.05
degrees_of_freedom = n - 1
t_score = qt(p=alpha/2, df=degrees_of_freedom,lower.tail=F)
margin_err_mae <- t_score * se_mae
margin_err_rmse <- t_score * se_rmse

low_bnd_mae <- mean_mae - margin_err_mae
up_bnd_mae <- mean_mae + margin_err_mae
print(c(low_bnd_mae,up_bnd_mae))
low_bnd_rmse <- mean_rmse - margin_err_rmse
up_bnd_rmse <- mean_rmse + margin_err_rmse
print(c(low_bnd_rmse,up_bnd_rmse))



