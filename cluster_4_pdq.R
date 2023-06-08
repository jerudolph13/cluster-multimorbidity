
###################################################################################################################
#
# Project: Multimorbidity clusters in ALIVE
#
# Purpose: Run probabilistic clustering
#
# Author: Jacqueline Rudolph
#
# Last Update: 06 June 2023
#
###################################################################################################################


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
  # Unlike other functions, PDQ takes binary variables as numeric
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

write_csv(cv.res, "../results/cv-metrics_pdq.csv")


# Run algorithm -----------------------------------------------------------

k = 2
set.seed(k)

# Probabilistic clustering
pdq <- PDQ(x=dat, k=k, ini="random", dist="gower", bin=1:ncol(dat))

# Obtain cluster assignment for each observation
clusters <- pdq$label

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
write_csv(dat3, file=paste0("../results/clusters_pdq", k, ".csv"))


# Examine cluster probabilities -------------------------------------------

all.prob <- pdq$probability
prob <- rep(NA, length(clusters))

# Examine distribution of probabilities by assigned cluster
for (i in 1:length(clusters)) {
  prob[i] <- all.prob[i, clusters[i]]
}

prob2 <- bind_cols(cluster=clusters, prob=prob) %>% 
  group_by(cluster) %>% 
  summarize(min = min(prob),
            p25 = quantile(prob, p=0.25),
            med = median(prob),
            mean = mean(prob),
            p75 = quantile(prob, p=0.75),
            max = max(prob))
