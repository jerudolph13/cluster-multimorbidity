
###########################################################################
#
# Project: Multimorbidity clusters in ALIVE
#
# Purpose: Run ensemble clustering
#
# Author: Jacqueline Rudolph, Bryan Lau
#
# Last Update: 16 Feb 2023
#
###########################################################################


packages <- c("tidyverse","data.table", "zoo", "Hmisc", "cluster", "fpc")
for (package in packages) {
  library(package, character.only=T)
}

source("./ensemble_util.R")
source("./cv_ensemble.R")


# Read in data ------------------------------------------------------------

# Read in data
dat.full <- read_csv(file="../data/alive_cluster_setup_age50.csv")

# Define variables to use when clustering
dat <- dat.full %>%
  select(hypertension, diabetes, lung_disease, renal_disease, liver_disease, 
         depression, cancer, hiv) %>% 
  mutate(across(everything(), ~ as.factor(.x)))

# Clustering results
hier2 <- read_csv(file="../results/clusters_hier2.csv") %>% 
  rename(hier2 = cluster)
hier3 <- read_csv(file="../results/clusters_hier3.csv") %>% 
  rename(hier3 = cluster)
pam2 <- read_csv(file="../results/clusters_pam2.csv") %>% 
  rename(pam2 = cluster)
pam3 <- read_csv(file="../results/clusters_pam3.csv") %>% 
  rename(pam3 = cluster)
pam5 <- read_csv(file="../results/clusters_pam5.csv") %>% 
  rename(pam5 = cluster)
pdq2 <- read_csv(file="../results/clusters_pdq2.csv") %>% 
  rename(pdq2 = cluster)
pdq5 <- read_csv(file="../results/clusters_pdq5.csv") %>% 
  rename(pdq5 = cluster)

cluster.mem <- select(hier2, id, hier2) %>% 
  left_join(select(hier3, id, hier3), by="id") %>% 
  left_join(select(pam2, id, pam2), by="id") %>% 
  left_join(select(pam3, id, pam3), by="id") %>% 
  left_join(select(pam5, id, pam5), by="id") %>% 
  left_join(select(pdq2, id, pdq2), by="id") %>% 
  left_join(select(pdq5, id, pdq5), by="id")


# Create list of hyperparameters ------------------------------------------

range.a1 <- seq(0.1, 1.0, 0.1)
range.a2 <- seq(0.1, 0.9, 0.1)

a1 <- rep(range.a1, length(range.a2)) %>% 
  sort(decreasing=T)
a2 <- rep(seq(0.1, 0.9, 0.1), length(range.a1))


# Determine optimal hyperparamters ----------------------------------------

# Compute cluster quality indices using cross-validation
cv.res <- mapply(x=a1, y=a2, 
              function(x, y){uml_cv(alpha1=x, alpha2=y, obs=dat, clust=select(cluster.mem, -id), nsplit=5)},
              SIMPLIFY=F)
cv.res <- do.call(rbind, cv.res)

write_csv(cv.res, "../results/cv-metrics_ensemble.csv")

# Rank hyperparameters
ps <- arrange(cv.res, desc(mean_ps))
ch <- arrange(cv.res, desc(mean_ch))
si <- arrange(cv.res, desc(mean_si))


# Run algorithm -----------------------------------------------------------

# Transform cluster membership matrix to binary indicator matrix
res1 <- cl.sim(tmp=select(cluster.mem, -id), bin.ind=T, wt=NULL)

# Merge at alpha1=0.8
set.seed(123)
res2 <- merge.clust(cls.res=res1, alpha=0.7)
  # Calculate membership similarity
  member <- memb.sim(zdat=res2)
  # Assign cluster based on membership similarity (alpha2=0.1)
  new.clust <- assign.clust(zdat=member, alpha2=0.1, wt=rep(1, nrow(member)))
    cluster <- new.clust[[1]] %>% 
      mutate(cluster = factor(cluster, levels=c(1,2,4,5), labels=c(1,2,3,4)))
    table(cluster)

  
# Merge and output --------------------------------------------------------

cluster.mem2 <- bind_cols(cluster.mem, ensemble=cluster$cluster)

dat2 <- left_join(dat.full, cluster.mem2, by="id")

write_csv(dat2, file="../results/res_cluster_all.csv")


