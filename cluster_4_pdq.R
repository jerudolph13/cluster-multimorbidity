
###################################################################################################################
#
# Project: Multimorbidity clusters in ALIVE
#
# Purpose: Run hierarchical clustering
#
# Author: Jacqueline Rudolph
#
# Last Update: 07 Feb 2023
#
###################################################################################################################

# Citation for PDQ method: https://www.sciencedirect.com/science/article/pii/S1568494622007530

packages <- c("tidyverse", "clue", "fpc", "zoo", "cluster", "FPDclustering")
for (package in packages) {
  library(package, character.only=T)
}

source("./cv_uml.R")


# Set up ------------------------------------------------------------------

set.seed(123)

# Read in data
dat.full <- read_csv(file="../data/alive_cluster_setup_age50.csv")

# Define variables to use when clustering
  # Unlike other functions, PDQ seems to take binary variables as numeric
dat <- dat.full %>%
  select(hypertension, diabetes, lung_disease, renal_disease, liver_disease, 
         depression, cancer, hiv) 


# Optimal number of clusters (k) ------------------------------------------

# Indices in full data
  # Unlike other functions, PDQ will not accept k=1
full.res <- lapply(2:10, function(x){uml_cv(k=x, data=dat, algorithm="pdq", nsplit=1)})
full.res <- do.call(rbind, full.res)

write_csv(full.res, "../results/full-metrics_pdq.csv")

# Cross-validation
cv.res <- lapply(2:10, function(x){uml_cv(k=x, data=dat, algorithm="pdq", nsplit=5)})
cv.res <- do.call(rbind, cv.res)

write_csv(cv.res, "../results/cv-metrics_pam.csv")


# Run algorithm -----------------------------------------------------------

k = 5
set.seed(k)

# Create dissimilarity matrix
diss.matrix <- daisy(dat, metric="gower")

# Hierarchical clustering using Complete Linkage
pam <- pam(diss.matrix, k=k, diss=T)

# Obtain cluster assignment for each observation
clusters <- pam$clustering

# Merge clusters back to data
dat2 <- bind_cols(dat.full, cluster=clusters) 

# Rank clusters by density of chronic conditions (from fewest conditions to most conditions)
ranks <- dat2 %>% 
  group_by(cluster) %>% 
  summarize(avg_hyper = mean(hypertension),
            avg_diab = mean(diabetes),
            avg_lung = mean(lung_disease),
            avg_renal = mean(renal_disease),
            avg_cirr = mean(liver_disease),
            avg_can = mean(cancer),
            avg_dep = mean(depression),
            avg_hiv = mean(hiv))

rowmeans <- rowMeans(select(ranks, -cluster))
ranks <- bind_cols(ranks, rowmeans=rowmeans) %>% 
  arrange(rowmeans) %>% 
  mutate(ranked_cluster=c(1:k))

dat3 <- merge(dat2, select(ranks, cluster, ranked_cluster), by="cluster") %>% 
  select(id, cluster, ranked_cluster) 

# Output assigned clusters for ensemble method
write_csv(dat3, file=paste0("../results/clusters_pam", k, ".csv"))




