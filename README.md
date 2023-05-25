# cluster-multimorbidity
This repository contains the R code for the in-progress manuscript "Characterizing multimorbidity in ALIVE: Comparing single and ensemble clustering methods." There are 3 utility functions:
  - cv_uml.R - runs crossvalidation for unsupervised machine learning algorithms across values of k
  - cv_ensemble.R - runs crossvalidation for the ensemble algorithm across values of alpha_1 and alpha_2
  - ensemble_util.R - runs the ensemble algorithm

We also include our programs setting up the data and implementing these functions
  - cluster_1_data.R - performs basic data management
  - cluster_2_impute.R - imputes missingness using random forest imputation
  - cluster_3_setup.R - sets up the imputed data for analysis
  - cluster_4_hier.R - runs hierarchical clustering
  - cluster_4_pam.R - runs PAM
  - cluster_4_pdq.R - runs probabilistic clustering
  - cluster_5_ensemble.R - uses the results from the single algorithms in the ensemble algorithm

