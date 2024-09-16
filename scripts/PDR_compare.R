###################
# Description: 
# This file compares the ability of clinician annotated Posterior Dominant Rhythm
# (PDR) vs automatically detected PDR to predict the age of individuals from our
# control cohort. We fit GAMs to predict age of an individual based on their PDR.
# This script reproduces results seen in Figure 2 of Galer et al.

###################


# Set working directory
setwd("/Volumes/helbig_lab/Users/galerp/EEG/manuscript/github/")

# Loads primary functions
source("main_R_functions.R")

# PDR information from clinicians (clin_dom_freq) and the various detected PDR
# via our automated method across iterations (frequency)
pdr_comb <- read_csv("data/pdr_controls_auto_vs_clin.csv")


##########
# Use filters to find median automatically detected PDR
##########

control_pdr_auto <- pdr_comb %>% 
  filter(frequency<=14) %>% 
  filter(frequency>=3) %>% 
  filter(prominence>5.5) %>%
  filter(power>6) %>% 
  filter(!(frequency<5 & power<13)) %>%
  filter(!(frequency>5 & width<8)) %>%
  filter(!(frequency>=9 & prominence<5.5)) %>%
  filter(!(frequency>=10 & prominence<7.25)) %>%
  filter(!(frequency>=11 & prominence<7.6)) %>%
  filter(width_height>1) %>% 
  filter(power<21) %>%
  filter(width<20) %>% 
  group_by(patient, age) %>%
  dplyr::summarise(frequency = median(frequency)) %>%
  rename(auto_freq = frequency)



#######################
# Compare PDRs
#######################
library(ggpubr)
pdr_med <- control_pdr_auto %>% 
  select(patient, auto_freq) %>% 
  distinct() %>% 
  left_join(pdr_comb %>% 
              select(patient, clin_dom_freq,age) %>% 
              distinct())%>% 
  filter(!is.na(age)) %>% 
  mutate(peak_dif_abs = abs(auto_freq - clin_dom_freq))%>% 
  mutate(peak_dif = auto_freq - clin_dom_freq)%>% 
  filter(!is.na(peak_dif)) %>% 
  rowwise() %>% 
  mutate(avg = mean(c(clin_dom_freq,auto_freq)))


# Directly compare PDRs and test Pearson correlation
ggplot(pdr_med, aes(x = clin_dom_freq, y = auto_freq)) +
  geom_point(size=2,fill = "aquamarine2", color="royalblue4", alpha = 0.4) +
  # geom_smooth(method='lm',aes(color = gene))+
  geom_abline(intercept = 0,linetype="dashed", linewidth=0.5) +
  coord_cartesian(xlim=c(3.1, 13),ylim=c(3.1, 13))+
  stat_cor(label.y.npc="top", label.x.npc = "left", method="pearson", size = 2.5)+
  ggtitle("Clinician vs Calculated EEG") +
  xlab("Clinician Annotation") +
  ylab("Calculated EEG")+
  theme_classic()+
  theme(legend.position = "none") +  
  theme_classic()+
  scale_x_continuous(breaks = seq(3, 13, by = 2)) +
  scale_y_continuous(breaks = seq(3, 13, by = 2)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))

# ggsave("PDR_manual_vs_clinician_cor.png",width = 6.5, height = 6, dpi = 1000)

###########
# Turn into long format
###########

pdr_long <- pdr_med %>%
  pivot_longer(
    cols = c(clin_dom_freq, auto_freq),  # Specify columns to pivot into longer format
    names_to = "Method",           # This will be the new column for the variable names (dom_freq, foof_peak)
    values_to = "PDR"           
  )
pdr_long$Method[pdr_long$Method=="clin_dom_freq"] = "Clinician"
pdr_long$Method[pdr_long$Method=="auto_freq"] = "Automated"



#########################
# Train Models
#########################


####################
# Model permutations
####################

pdr_auto <- pdr_long %>% 
  filter(Method == "Automated")

pdr_clin <- pdr_long %>% 
  filter(Method == "Clinician")


# Initialize the results dataframe
results_clin_df <- data.frame()
results_auto_df <- data.frame()

# How large training set it
ntrain = round(length(unique(pdr_long$patient))*0.8)  


# Perform 1000 permutations with 80-20 training-testing split
for (i in 1:1000) {
  
  train_pats = sample(unique(pdr_long$patient), ntrain)
  # Separate training and testing data
  pdr_auto_train <- pdr_auto %>% 
    filter(patient %in% train_pats)
  pdr_auto_test <- pdr_auto %>% 
    filter(patient %nin% train_pats)
  
  pdr_clin_train <- pdr_clin %>% 
    filter(patient %in% train_pats)
  pdr_clin_test <- pdr_clin %>% 
    filter(patient %nin% train_pats)
  
  
  
  gam_model <- gam(age ~ s(PDR), data = pdr_clin_train, family = poisson(link = "log"))
  pdr_clin_test$PDR_pred <- predict(gam_model, newdata = pdr_clin_test, type = "response")
  # Calculate metrics
  abs_errors = abs(pdr_clin_test$PDR - pdr_clin_test$PDR_pred)
  mse <- mean(abs_errors)
  rmse <- sqrt(mean((pdr_clin_test$PDR - pdr_clin_test$PDR_pred)^2))
  r_squared <- cor(pdr_clin_test$PDR,  pdr_clin_test$PDR_pred)^2
  ci <- quantile(abs_errors, probs = c(0.025, 0.975))
  # Store the results
  results_clin_df <- rbind(results_clin_df, data.frame(Iteration = i, MSE = mse, RMSE = rmse, R_Squared = r_squared, ci_lower = ci[1], ci_upper = ci[2], ci_width = (ci[2]-ci[1])))
  
  
  gam_model <- gam(age ~ s(PDR), data = pdr_auto_train, family = Gamma(link = "log"))
  pdr_auto_test$PDR_pred <- predict(gam_model, newdata = pdr_auto_test, type = "response")
  # Calculate metrics
  abs_errors = abs(pdr_auto_test$PDR - pdr_auto_test$PDR_pred)
  mse <- mean(abs_errors)
  rmse <- sqrt(mean((pdr_auto_test$PDR - pdr_auto_test$PDR_pred)^2))
  r_squared <- cor(pdr_auto_test$PDR,  pdr_auto_test$PDR_pred)^2
  ci <- quantile(abs_errors, probs = c(0.025, 0.975))
  # Store the results
  results_auto_df <- rbind(results_auto_df, data.frame(Iteration = i, MSE = mse, RMSE = rmse, R_Squared = r_squared, ci_lower = ci[1], ci_upper = ci[2], ci_width = (ci[2]-ci[1])))
  
}

results_comp <- results_clin_df %>% 
  mutate(model = "Clinician") %>% 
  rbind(results_auto_df %>% 
          mutate(model = "Automated"))


results_dif <- results_clin_df %>% 
  rename(clin_MSE = MSE, clin_RMSE = RMSE, clin_R_Squared = R_Squared, clin_ci_lower = ci_lower, clin_ci_upper = ci_upper, clin_ci_width = ci_width) %>% 
  left_join(results_auto_df %>% 
              rename(auto_MSE = MSE, auto_RMSE = RMSE, auto_R_Squared = R_Squared, auto_ci_lower = ci_lower, auto_ci_upper = ci_upper, auto_ci_width = ci_width)) %>% 
  mutate(dif_MSE = clin_MSE - auto_MSE, dif_RMSE = clin_RMSE - auto_RMSE )



######
# Find iteration with median RMSE difference
######
rmse_dif_med <- results_dif %>% 
  mutate(dif_RMSE_med = abs(dif_RMSE-median(results_dif$dif_RMSE))) %>% 
  arrange(dif_RMSE_med) %>% 
  select(Iteration, dif_RMSE_med)

iter_rmse_dif_med <- rmse_dif_med$Iteration[1]

####################
# Plot results
####################


# Mean squared Error
ggplot(results_dif, aes(auto_MSE, clin_MSE))+
  geom_point(color = "purple4", alpha = 0.55, size = 2)+
  geom_abline(intercept = 0,linetype="dashed", linewidth=0.75) +
  coord_cartesian(xlim=c(1.9, 3.3),ylim=c(1.9, 3.3))+
  xlab("Automated Model MAE")+
  ylab("Clinician Model MAE")+
  # theme_minimal() +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))

# Root mean squared error
ggplot(results_dif, aes(auto_RMSE, clin_RMSE))+
  geom_point(color = "purple4", alpha = 0.55, size = 2)+
  geom_abline(intercept = 0,linetype="dashed", linewidth=0.75) +
  coord_cartesian(xlim=c(2.25, 3.6),ylim=c(2.25, 3.6))+
  xlab("Automated Model RMSE")+
  ylab("Clinician Model RMSE")+
  # theme_minimal() +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))



# Transform for next figures
results_long <- results_dif %>% 
  select(Iteration, clin_MSE, clin_RMSE) %>% 
  rename(MSE = clin_MSE, RMSE = clin_RMSE) %>% 
  mutate(Model = "Clinician") %>% 
  rbind(results_dif %>% 
          select(Iteration, auto_MSE, auto_RMSE) %>% 
          rename(MSE = auto_MSE, RMSE = auto_RMSE) %>% 
          mutate(Model = "Automated"))

# Figure seen in paper (Figure 2C Inset)
# ROOT MEAN SQUARED ERROR COMPARISON FIGURE
results_long$Model <- factor(results_long$Model, levels = c("Clinician", "Automated"))


RMSE_med <- results_long %>% 
  filter(Iteration == iter_rmse_dif_med)


ggplot(results_long, aes(x = Model, y = RMSE, group = Iteration)) +
  geom_line(aes(group = Iteration), color = "grey68", alpha = 0.1) +
  # geom_point(aes(color = Model), size = 2, alpha =0.1) +
  geom_jitter(aes(color = Model), position = position_dodge(width = 0.3), alpha = 0.2) +
  
  geom_line(data = results_long %>% filter(Iteration == iter_rmse_dif_med),
            aes(x = Model, y = RMSE),
            color = "red3", size = 1.5, alpha = 0.6)+
  ylab("RMSE")+
  theme_minimal() +
  scale_color_manual(values = c("Automated" = "#F8726A", "Clinician" = "turquoise3"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))


ggsave("PDR_manual_vs_clin_RMSE_dot_line.png", width = 6, height = 5.5, dpi = 1000)



