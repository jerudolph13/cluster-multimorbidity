
###################################################################################################################
#
# Project: Multimorbidity clusters in ALIVE
#
# Purpose: Run hierarchical clustering
#
# Author: Jacqueline Rudolph
#
# Last Update: 13 Feb 2023
#
###################################################################################################################


packages <- c("tidyverse", "clue", "fpc", "zoo", "cluster", "NbClust")
for (package in packages) {
  library(package, character.only=T)
}

source("./cv_uml.R")


# Set up ------------------------------------------------------------------

# Read in data
dat.full <- read_csv(file="../data/alive_cluster_setup_age50.csv")

# Define variables to use when clustering
dat <- dat.full %>%
  select(hypertension, diabetes, lung_disease, renal_disease, liver_disease, 
         depression, cancer, hiv) %>% 
  mutate(across(everything(), ~ as.factor(.x)))


# Optimal number of clusters (k) ------------------------------------------

# Indices in full data
full.res <- lapply(1:10, function(x){uml_cv(k=x, data=dat, algorithm="hierarchical", nsplit=1)})
full.res <- do.call(rbind, full.res)
  # k = 2 wins based on CH and SI

write_csv(full.res, "../results/full-metrics_hier.csv")

# Cross-validation
cv.res <- lapply(1:10, function(x){uml_cv(k=x, data=dat, algorithm="hierarchical", nsplit=5)})
cv.res <- do.call(rbind, cv.res)
  # k = 2 wins based on PS (but value below 0.8)
  # k = 2 wins based on CH
  # k = 3 wins based on SI

write_csv(cv.res, "../results/cv-metrics_hier.csv")


# Run algorithm -----------------------------------------------------------

k = 3
set.seed(k)

# Create dissimilarity matrix
diss.matrix <- daisy(dat, metric="gower")

# Hierarchical clustering using Complete Linkage
hier <- hclust(diss.matrix, method="complete")

# Obtain cluster assignment for each observation
clusters <- tibble(cutree(hier, k=k))
names(clusters) <- "cluster"

# Merge clusters back to data
dat2 <- bind_cols(dat.full, clusters) 

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
write_csv(dat3, file=paste0("../results/clusters_hier", k, ".csv"))




