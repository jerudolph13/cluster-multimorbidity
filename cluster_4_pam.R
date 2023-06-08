
###################################################################################################################
#
# Project: Multimorbidity clusters in ALIVE
#
# Purpose: Run partition around medoids
#
# Author: Jacqueline Rudolph
#
# Last Update: 07 Feb 2023
#
###################################################################################################################


packages <- c("tidyverse", "clue", "fpc", "zoo", "cluster")
for (package in packages) {
  library(package, character.only=T)
}

source("./cv_uml.R")


# Set up ------------------------------------------------------------------

set.seed(123)

# Read in data
dat.full <- read_csv(file="../data/alive_cluster_setup_age50.csv")

# Define variables to use when clustering
dat <- dat.full %>%
  select(hypertension, diabetes, lung_disease, renal_disease, liver_disease, 
         depression, cancer, hiv) %>% 
  mutate(across(everything(), ~ as.factor(.x)))


# Optimal number of clusters (k) ------------------------------------------

# Compute cluster quality indices in full data
full.res <- lapply(1:10, function(x){uml_cv(k=x, data=dat, algorithm="pam", nsplit=1)})
full.res <- do.call(rbind, full.res)

write_csv(full.res, "../results/full-metrics_pam.csv")

# Compute cluster quality indices using cross-validation
cv.res <- lapply(1:10, function(x){uml_cv(k=x, data=dat, algorithm="pam", nsplit=5)})
cv.res <- do.call(rbind, cv.res)

write_csv(cv.res, "../results/cv-metrics_pam.csv")


# Run algorithm -----------------------------------------------------------

k = 5
set.seed(k)

# Create dissimilarity matrix
diss.matrix <- daisy(dat, metric="gower")

# Partition around medoids
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




