
###########################################################################
#
# Project: Multimorbidity clusters
#
# Purpose: Impute data
#
# Author: Jacqueline Rudolph
#
# Last Update: 01 Feb 2023
#
###########################################################################


packages <- c("tidyverse", "lubridate", "tidyselect", "missForest", "doParallel")
for (package in packages) {
  library(package, character.only=T)
}


# Read in data ------------------------------------------------------------

dat <- read_csv("../data/alive_cluster_long.csv") %>% 
  mutate(across(!c(id, visdate, bpsys, bpdias, hba1c, fev_fvc, gfr, 
                   uprt_crt, fbscan, cesdtot, age, n_vis), ~ factor(.x)))


# Manage data -------------------------------------------------------------

# Pull out baseline variables
base <- select(dat, id, m0f1, black, enroll_period, beduchs) %>% 
  filter(!duplicated(id))

# Make data set wide
wide <- dat %>% 
  mutate(yr = year(visdate)) %>% 
  select(-c(visdate, m0f1, black, enroll_period, beduchs)) %>% 
  pivot_wider(names_from=n_vis, values_from=-c(id, n_vis))

# Merge baseline data back
wide2 <- as.data.frame(left_join(base, wide, by="id"))


# Impute data -------------------------------------------------------------

set.seed(123)
registerDoParallel(cores=7) # Run in parallel for faster computation
impute <- missForest(select(wide2, -id), verbose=T, ntree=100, parallelize="forests")
wide.impute <- impute$ximp
error <- impute$OOBerror

wide.impute2 <- bind_cols(id=wide2$id, wide.impute)


# Make data long ----------------------------------------------------------

long <- wide.impute2 %>% 
  pivot_longer(cols=!c("id", "m0f1", "black", "enroll_period", "beduchs"), 
               names_to = c(".value", "n_vis"),
               names_pattern = "(.*)_(.*)") %>% 
  mutate(n_vis=as.numeric(n_vis))

# Keep imputed values from observed visits
long2 <-  dat %>% 
  select(id, visdate, n_vis) %>% 
  left_join(long, by=c("id", "n_vis"))


# Output data -------------------------------------------------------------

write_csv(long2, "../data/alive_cluster_imputed.csv")
