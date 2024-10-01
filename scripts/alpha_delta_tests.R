# These are packages only required for this R script
library(gamm4) #to create GAMM models
library(pbkrtest) #to test fit difference between GAMMs
library(FSA) #for Dunn test

###################
# Description: 
# This file plots alpha-delta ratio between genes and variants and performs the 
# relevant statistical analyses. This code should reproduce results from Figure 2.

###################

# Set working directory
setwd("/Volumes/helbig_lab/Users/galerp/EEG/manuscript/scalp_EEG_DEE/")

# Loads primary functions
source("scripts/main_R_functions.R")

# Spectral information for every electrode
psd_bp_rel <- read_csv("data/psd_bp_gene_controls.csv")

# STXBP1 variant information
stx_vars <- read_csv("data/stxbp1_var_type.csv")

####################
# Median Across Scalp
####################

psd_bp_rel_pat <- psd_bp_rel %>% 
  group_by(patient, age, age_group, gene, bandpower) %>% 
  dplyr::summarise(Power = median(power)) %>% 
  filter(age_group !="Neonates")

####################
# Alpha-delta ratio 
####################

psd_bp_ad_rat <- psd_bp_rel_pat %>% 
  filter(bandpower == "alpha_delta")

####################
# Alpha-delta ratio by STXBP1 Variant
####################
psd_bp_ad_vars <- psd_bp_ad_rat %>% 
  left_join(stx_vars)%>%
  mutate(var_group = case_when(gene=="Controls"~"Controls",
                               (is.na(var_group)&gene=="SCN1A")~"NA_SCN1A",
                               (is.na(var_group)&gene=="SYNGAP1")~"NA_SYNGAP1",
                               TRUE~paste(gene,var_group,sep="_")))


ggplot(psd_bp_ad_vars %>%
         filter(gene %in% c("STXBP1", "Controls")) , 
       aes(x = age, y = Power, fill = var_group, color = var_group)) +
  geom_smooth(data = . %>% filter(var_group != "NA_STXBP1")) +  # Plots all data points including "Controls"
  geom_point(data = . %>% filter(gene != "Controls")) +  # Excludes "Controls" for plotting points
  # geom_point(data = . %>% filter(gene == "Controls"), alpha = 0.1) +  # Excludes "Controls" for plotting points
  
  coord_cartesian(xlim=c(0, 18),ylim=c(-0.04, 0.58))+
  scale_fill_manual(breaks = c("STXBP1_PTV","STXBP1_missense","Controls"),
                    values = c("#C77CFF", "#00BFC4","#F8766D"))+
  scale_color_manual(breaks = c("STXBP1_PTV", "STXBP1_missense","Controls"),
                     values = c("#C77CFF", "#00BFC4","#F8766D"))+
  xlab("Age")+
  ylab("Alpha Delta Power")+
  ggtitle("Alpha Delta Power in STXBP1 Variants")+
  theme_classic()+ 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))

# ggsave("alpha_delta_stx_var_wide.png",  width = 9, height = 4.5, dpi=1000)


######
# Now test if significant with GAMM test
######

psd_stx <- psd_bp_ad_vars %>%
  filter(gene =="STXBP1")

psd_stx$var_group <- as.factor(psd_stx$var_group)
psd_stx$patient <- as.factor(psd_stx$patient)


# Fit models with GAMM
model_0 <- gamm4(Power ~ s(age), data = psd_stx, random = ~(1|patient)) #"Null" model
model_1 <- gamm4(Power ~ s(age) + var_group, data = psd_stx, random = ~(1|patient)) #"Variant" model

# Returns whether model with the variant group has a significantly better fit
KRmodcomp(model_1$mer,model_0$mer)$test$p.value[1]


####################
# Now alpha-delta difference between genes
## By age group
####################

psd_bp_ad_rat$gene <- factor(psd_bp_ad_rat$gene, 
                             levels=c("Controls", "SCN1A", "SYNGAP1","STXBP1"))
# Colors are kept consistent across the paper
my_colors <- c(Controls = "#F8766D", SCN1A= "#7CAE00", STXBP1 = "#C77CFF", SYNGAP1 = "#00BFC4")

desired_order <- c("Infant","Toddler", "Child","Adult")


psd_bp_ad_rat$age_group <- factor(psd_bp_ad_rat$age_group, levels = desired_order)


ggplot(psd_bp_ad_rat %>% filter(age<=21), aes(age_group, Power))+
  geom_boxplot(aes(fill = gene,color = gene),alpha=0.25,outliers = T)+
  # To aid in visualization, removing outliers 
  # NOTE this function does not remove from calculation of boxplot, just removes from visualization
  # Do not use ylim(c()) as this removes from boxplot calculation
  coord_cartesian(ylim=c(0, 1.1)) +
  scale_y_continuous(breaks = round(seq(0,1.12, by = 0.25),2))+
  ggtitle("Alpha-Delta Ratio")+
  ylab("Alpha-Delta Ratio")+
  xlab("Age Group")+
  theme_classic() + 
  scale_color_manual(values = my_colors)+
  theme(text = element_text(size = 14),axis.text = element_text(size = 14))
 
# ggsave("alpha_delta_all_age_boxplot_21y.png", width=9.5, height = 6, dpi = 1000)

#####
# Test statistical difference between genes by age group
#####
# Computing Dunn test with Benjamini-Hochberg correction and computing Cohen's d

# For Genes
med_dunn_test <- as.data.frame(matrix(nrow= 0, ncol = 7))
names(med_dunn_test) <- c("age_group", "Comparison", "feature", "Z","p_unadj","p_adj")


psd_bp_ad_rat_dunn <- psd_bp_ad_rat %>% 
  #Triple checking
  filter(age_group !="Neonates") %>% 
  filter(age<=21)



for (a in unique(psd_bp_ad_rat_dunn$age_group)) {
  psd_sub <- psd_bp_ad_rat_dunn %>%
    filter(age_group == a)
  
  # Extract group statistics for each gene group
  ## This is to calculate Cohen's d later
  group_stats <- psd_sub %>%
    group_by(gene) %>%
    dplyr::summarize(
      mean = mean(Power),
      sd = sd(Power),
      n = n()
    )
  
  d_test <- dunnTest(Power ~ gene, data = psd_sub, method = "bh")
  d_test_df <- d_test$res %>% as.data.frame()
  d_test_df$age_group <- a
  
  # Calculate Cohen's d for each gene comparison in d_test_df
  for (i in 1:nrow(d_test_df)) {
    # Extract group names from comparison ("GeneA - GeneB")
    comparison <- strsplit(d_test_df$Comparison[i], " - ")[[1]]
    gene_A <- comparison[1]
    gene_B <- comparison[2]
    
    # Get the means, SDs, and sample sizes for both gene groups
    mean_A <- group_stats %>% filter(gene == gene_A) %>% pull(mean)
    mean_B <- group_stats %>% filter(gene == gene_B) %>% pull(mean)
    sd_A <- group_stats %>% filter(gene == gene_A) %>% pull(sd)
    sd_B <- group_stats %>% filter(gene == gene_B) %>% pull(sd)
    n_A <- group_stats %>% filter(gene == gene_A) %>% pull(n)
    n_B <- group_stats %>% filter(gene == gene_B) %>% pull(n)
    
    #  Pooled standard deviation
    pooled_sd <- sqrt(((n_A - 1) * sd_A^2 + (n_B - 1) * sd_B^2) / (n_A + n_B - 2))
    
    #  Cohen's d
    cohens_d <- (mean_A - mean_B) / pooled_sd
    
    # Add Cohen's d to the current row 
    d_test_df$Cohen_d[i] <- cohens_d
  }
  
  # Bind Cohen's d to med_dunn_test
  med_dunn_test <- med_dunn_test %>% rbind(d_test_df)
}

print(med_dunn_test)


