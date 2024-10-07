# Scalp EEG Analysis in Childhood Genetic Epilepsies



This repository is composed of a series of scripts which analyzes the spectral features from 2036 clinical scalp EEGs from 1957 EEGs from 1585 individuals ages 0.03-38.68 years (median 7.23 years, IQR 2.61-12.37 years). Most analyses cover three genetic epilepsies and a large control cohort: _STXBP1_ (95 EEGs from 20 individuals, ages 0.16-17.77 years), _SCN1A_ (154 EEGs from 68 individuals, ages 0.30-24.62 years), _SYNGAP1_ (46 EEGs from 21 individuals, ages 1.19-27.47 years), and controls (847 EEGs from 806 individuals, ages 0.03-38.68 years). There are also two additionial cohorts: individuals with a scalp EEG and a seizure frequency annotation collected close in time (440 EEGs from 354 individuals, ages 1.59-20.03 years) and individuals with a scalp EEG and gross motor function assessment (GMFM) collected close in time (400 EEGs from 340 individual, ages 0.64 to 22.18 years).
These analyses should replicated the primary findings in the manuscript by Galer et al. "Quantitative EEG Spectral Features Differentiate Genetic Epilepsies and Predict Neurologic Outcomes".
EEG cleaning and extraction of spectral features and posterior dominant rhythm were computed in Python. Downstream analyses including machine learning models and statistical analyses were computed in R.
An example EEG deidentified file from a control patient in our cohort can be retrieved here: https://upenn.box.com/v/galerqeegcntrleg. For the purposes of deidentifications, however, the original edf file was converted to a fif file which has removed some information which was used in our pipeline for cleaning (e.g., clinician annotations).


## Scripts:

* [Main R Functions](https://github.com/galerp/scalp_EEG_DEE/tree/main/scripts/main_R_functions.R) - Main R functions used throughout many of the scripts.
  
* [PDR Analysis](https://github.com/galerp/scalp_EEG_DEE/tree/main/scripts/PDR_compare.R) - extracts posterior dominant rhythms (PDR) from controls and compares it to clinician annotated PDR. Trains models using each source of PDR to predict age of the individual and compares the accuracy.
  
* [Alpha-Delta Ratio Analyses](https://github.com/galerp/scalp_EEG_DEE/tree/main/scripts/alpha_delta_tests.R)  - Compares the alpha-delta ratio across different gene groups and compared to controls.
  
* [Gene Prediction](https://github.com/galerp/scalp_EEG_DEE/tree/main/scripts/gene_prediction.R) - Trains and tests random forest models using spatial and non-spatial spectral features from EEG to differentiate controls from specific gene populations. Also trains and tests and three-class model with all three gene groups.

* [Gene Prediction](https://github.com/galerp/scalp_EEG_DEE/tree/main/scripts/gene_prediction.R) - Trains and tests random forest models using spatial and non-spatial spectral features from EEG to differentiate controls from specific gene populations. Also trains and tests and three-class model with all three gene groups.

* [Seizure Frequency Prediction](https://github.com/galerp/scalp_EEG_DEE/tree/main/scripts/seiz_freq_pred.R) - Trains and tests random forest models using spectral features from EEG to predict seizure frequency of individuals.

* [GMFM Prediction](https://github.com/galerp/scalp_EEG_DEE/tree/main/scripts/GMFM_prediction.R) - Trains and tests random forest models using spectral features from EEG to predict gross motor function measure (GMFM) scores. Results are compared agains a null model trained on just age of the individual at the time of the GMFM.

* [EEG Cleaning and Spectral Extraction](https://github.com/galerp/scalp_EEG_DEE/tree/main/scripts/python_functions) - This file contains the main python functions used to clean and extract spectral features from the scalp EEGs.

* [EEG Cleaning and Spectral Extraction Walkthrough](https://github.com/galerp/scalp_EEG_DEE/blob/main/scripts/eeg_psd_wt.ipynb) -  This is a jupyter notebook which walks through the basic steps of cleaning an EEG file, removing artifact-heavy epochs, and extracting the power spectral density of each electrode. It requires example data that can be retrieved here: https://upenn.box.com/v/galerqeegcntrleg
  
## Files: ##

* [Controls' PDR](https://github.com/galerp/scalp_EEG_DEE/tree/main/data/pdr_controls_auto_vs_clin.csv)  - This file contains the output from our posterior dominant rhythm (PDR) detector across EEG epochs (prior to additionial filters) and the PDR annotated by clincians on the same respective EEG.

* [Main Cohort Spectral Features](https://github.com/galerp/scalp_EEG_DEE/tree/main/data/psd_bp_gene_controls.csv)  - This file contains the relevant spectral features for each electrode across the main cohorts, _STXBP1_, _SYNGAP1_, _SCN1A_, and Controls. Features are the median of each electrode across all available EEG epochs.

* [STXBP1 Variant Types](https://github.com/galerp/scalp_EEG_DEE/tree/main/data/stxbp1_var_type.csv)  - This file contains the broad category of variant type (protien truncating variant (PTV) vs missense variant) in individuals with a variant in _STXBP1_ in our cohort.

* [Seizure Frequency Spectral Features](https://github.com/galerp/scalp_EEG_DEE/tree/main/data/psd_bp_seiz_freq.csv)  - This file contains the relevant spectral features and respective seizure frequency for all individuals in the seizure frequency analyses.

* [GMFM Spectral Features](https://github.com/galerp/scalp_EEG_DEE/tree/main/data/psd_bp_GMFM.csv)  - This file contains the relevant spectral features and respective gross motor function measure (GMFM) score for all individuals in the GMFM analyses.


### R Requirements:
  These scripts use [R](https://www.r-project.org/) version 4.4.0 with the following packages:
1. [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html)
2. [hmisc](https://cran.r-project.org/web/packages/hmisc/index.html)
3. [reticulate](https://cran.r-project.org/web/packages/reticulate/index.html)
4. [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
5. [randomForest](https://cran.r-project.org/web/packages/randomForest/index.html)
6. [pROC](https://cran.r-project.org/web/packages/pROC/index.html)
7. [caret](https://cran.r-project.org/web/packages/caret/index.html)
8. [gamm4](https://cran.r-project.org/web/packages/gamm4/index.html)
9. [pbkrtest](https://cran.r-project.org/web/packages/pbkrtest/index.html)
10. [mgcv](https://cran.r-project.org/web/packages/mgcv/index.html)

### Python Dependencies:
These scripts use Python 3.9.13. The following third-party packages are required to run the Python scripts:
1.  [scipy](https://docs.scipy.org/doc/scipy/reference/signal.html](https://docs.scipy.org/doc/scipy/)
2.  [pandas](https://pandas.pydata.org/docs/)
3.  [mne](https://mne.tools/stable/index.html)
4.  [mne_icalabel](https://mne.tools/mne-icalabel/stable/index.html)

