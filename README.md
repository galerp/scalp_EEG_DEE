# Scalp EEG Analysis in Childhood Genetic Epilepsies


This repository is composed of a series of scripts which analyzes the spectral features from 2036 clinical scalp EEGs from 1676 individuals ages 0.03-38.71 years (median 7.38 years, IQR 2.59-12.49 years). Most analyses cover three genetic epilepsies and a large control cohort: _STXBP1_ (95 EEGs from 20 individuals, ages 0.16-17.77 years), _SCN1A_ (155 EEGs from 68 individuals, ages 0.30-24.62 years), _SYNGAP1_ (46 EEGs from 21 individuals, ages 1.19-27.47 years), and controls (2036 EEGs from 1676 individuals, ages 0.03-38.71 years). There are also two additionial cohorts: individuals with a scalp EEG and a seizure frequency annotation collected close in time (440 EEGs from 354 individuals, ages 1.59-20.03 years) and individuals with a scalp EEG and gross motor function assessment (GMFM) collected close in time (400 EEGs from 340 individual, ages 0.64 to 22.18 years).
These analyses should replicated the primary findings in the manuscript by Galer et al. "Quantitative EEG Spectral Features Differentiate Genetic Epilepsies and Predict Neurologic Outcomes".


## Scripts:

* [Binned Fishers test](https://github.com/galerp/Cube3/blob/main/scripts/fisher_dx_binned.R) - finds clinical associations with genes prior to diagnosis ("note elimination") in 3 month time bins.
* [Random Forest Models](https://github.com/galerp/Cube3/blob/main/scripts/rf_dx_model.R) - trains and tests Random Forest models to predict SCN1A, bootstrapping at every age interval.
* [HPO Propagation](https://github.com/galerp/Cube3/blob/main/additional_analyses/compose_prop.R)  - propagates a base HPO file, including all ancestors of each HPO in every time bin.


## Files: ##

[hpo def](https://github.com/galerp/Cube3/blob/main/Files/HPO_def_rl_2020-10-12_dl_2021-08-03.csv) - This file contains every HPO code along with its definition.
