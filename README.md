# Scalp EEG Analysis in Childhood Genetic Epilepsies


This repository is composed of a series of scripts which analyzes the spectral features from 2036 clinical scalp EEGs from 1676 individuals ages 0.03-38.71 years (median 7.38 years, IQR 2.59-12.49 years). Most analyses cover three genetic epilepsies and a large control cohort: _STXBP1_ (95 EEGs from 20 individuals, ages 0.16-17.77 years), _SCN1A_ (155 EEGs from 68 individuals, ages 0.30-24.62 years), _SYNGAP1_ (46 EEGs from 21 individuals, ages 1.19-27.47 years), and controls (2036 EEGs from 1676 individuals, ages 0.03-38.71 years). There are also two additionial cohorts: individuals with a scalp EEG and a seizure frequency annotation collected close in time (440 EEGs from 354 individuals, ages 1.59-20.03 years) and individuals with a scalp EEG and gross motor function assessment (GMFM) collected close in time (400 EEGs from 340 individual, ages 0.64 to 22.18 years).
These analyses should replicated the primary findings in the manuscript by Galer et al. "Quantitative EEG Spectral Features Differentiate Genetic Epilepsies and Predict Neurologic Outcomes".


## Scripts:

* [Main R Functions](https://github.com/galerp/scalp_EEG_DEE/tree/main/scripts/main_R_functions.R) - Main R functions used throughout many of the scripts.
  
* [PDR Analysis](https://github.com/galerp/scalp_EEG_DEE/tree/main/scripts/PDR_compare.R) - extracts posterior dominant rhythms (PDR) from controls and compares it to clinician annotated PDR. Trains models using each source of PDR to predict age of the individual and compares the accuracy.
  
* [Alpha-Delta Ratio Analyses](https://github.com/galerp/scalp_EEG_DEE/tree/main/scripts/alpha_delta_tests.R)  - Compares the alpha-delta ratio across different gene groups and compared to controls.
  
* [Gene Prediction](https://github.com/galerp/scalp_EEG_DEE/tree/main/scripts/gene_prediction.R))  - Trains and tests random forest models using spatial and non-spatial spectral features from EEG to differentiate controls from specific gene populations. Also trains and tests and three-class model with all three gene groups.

* [Gene Prediction](https://github.com/galerp/scalp_EEG_DEE/tree/main/scripts/gene_prediction.R))  - Trains and tests random forest models using spatial and non-spatial spectral features from EEG to differentiate controls from specific gene populations. Also trains and tests and three-class model with all three gene groups.

* [Seizure Frequency Prediction](https://github.com/galerp/scalp_EEG_DEE/tree/main/scripts/seiz_freq_pred.R)  - Trains and tests random forest models using spectral features from EEG to predict seizure frequency of individuals.

* [GMFM Prediction](https://github.com/galerp/scalp_EEG_DEE/tree/main/scripts/GMFM_prediction.R))  - Trains and tests random forest models using spectral features from EEG to predict gross motor function measure (GMFM) scores. Results are compared agains a null model trained on just age of the individual at the time of the GMFM.
  
## Files: ##

[hpo def](https://github.com/galerp/Cube3/blob/main/Files/HPO_def_rl_2020-10-12_dl_2021-08-03.csv) - This file contains every HPO code along with its definition.
