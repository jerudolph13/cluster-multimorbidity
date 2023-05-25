
###########################################################################
#
# Project: Multimorbidity clusters
#
# Purpose: Define chronic conditions & determine year of diagnosis
#
# Author: Jacqueline Rudolph
#
# Last Update: 01 Feb 2023
#
###########################################################################

packages <- c("tidyverse", "lubridate", "zoo")
for (package in packages) {
  library(package, character.only=T)
}


# Read in data ------------------------------------------------------------

dat <- read_csv("../data/alive_cluster_imputed.csv")

cancer <- read_csv("../data/alive_cancer.csv") %>% 
  select(id, dx_date) %>% 
  # Grab first cancers only
  filter(!duplicated(id))


# Create indicators of disease --------------------------------------------

dat <- left_join(dat, cancer, by="id")

dat2 <- dat %>% 
  # Diabetes: HbA1c > 6.5% or medication use
  mutate(hi_hba1c = as.numeric((hba1c>6.5)),
         
         # Obstructive lung disease: FEV/FVC <=0.70
         poor_lung_func = as.numeric((fev_fvc<=0.70)),
         
         # Kidney dysfunction: urine protein-creatinine ratio >200 or GFR <60
         renal_dysfunc = as.numeric((uprt_crt>200) | (gfr<60)),
         
         # Liver fibrosis: liver stiffness >=9.3
         fibrosis = as.numeric(fbscan>=9.3),
         
         # Liver cirrhosis: liver stiffness >=12.3
         cirrhosis = as.numeric(fbscan>=12.3),
         
         # Hypertension: Dias BP >=90, Sys BP >=140, or treatment in last 6 months
         hi_bp = as.numeric((bpdias>=90) | (bpsys>=140)),
         
         # Cancer: any diagnosis in yr or before
         cancer = ifelse(is.na(dx_date), 0, as.numeric(dx_date<=visdate)),
         
         # Depression: CESD>=23 or self-reported treatment in last 6 months
         cesd23 = as.numeric(cesdtot>=23))


# Calculate time between elevated measurements ----------------------------

hba1c <- dat2 %>% 
  filter(hi_hba1c==1) %>% 
  group_by(id) %>% 
  mutate(t_between = (lag(visdate) %--% visdate)/days(1),
         flag = t_between<365 & !is.na(t_between),
         diab_idx = cumsum(cumsum(flag))) %>% 
  filter(diab_idx==1) %>% 
  select(id, visdate, diab_idx) 

hbp <- dat2 %>% 
  filter(hi_bp==1) %>% 
  group_by(id) %>% 
  mutate(t_between = (lag(visdate) %--% visdate)/days(1),
         flag = t_between<365 & !is.na(t_between),
         hbp_idx = cumsum(cumsum(flag))) %>% 
  filter(hbp_idx==1) %>% 
  select(id, visdate, hbp_idx) 

lung_func <- dat2 %>% 
  filter(poor_lung_func==1) %>% 
  group_by(id) %>% 
  mutate(t_between = (lag(visdate) %--% visdate)/days(1),
         flag = t_between<365 & !is.na(t_between),
         lung_idx = cumsum(cumsum(flag))) %>% 
  filter(lung_idx==1) %>% 
  select(id, visdate, lung_idx) 

renal_func <- dat2 %>% 
  filter(renal_dysfunc==1) %>% 
  group_by(id) %>% 
  mutate(t_between = (lag(visdate) %--% visdate)/days(1),
         flag = t_between<365 & !is.na(t_between),
         renal_idx = cumsum(cumsum(flag))) %>% 
  filter(renal_idx==1) %>% 
  select(id, visdate, renal_idx) 

liver_fib <- dat2 %>% 
  filter(fibrosis==1) %>% 
  group_by(id) %>% 
  mutate(t_between = (lag(visdate) %--% visdate)/days(1),
         flag = t_between<365 & !is.na(t_between),
         liver_idx = cumsum(cumsum(flag))) %>% 
  filter(liver_idx==1) %>% 
  select(id, visdate, liver_idx) 

cesd <- dat2 %>% 
  filter(cesd23==1) %>% 
  group_by(id) %>% 
  mutate(t_between = (lag(visdate) %--% visdate)/days(1),
         flag = t_between<365 & !is.na(t_between),
         dep_idx = cumsum(cumsum(flag))) %>% 
  filter(dep_idx==1) %>% 
  select(id, visdate, dep_idx) 

dat3 <- dat2 %>% 
  left_join(hba1c, by=c("id", "visdate")) %>% 
  left_join(hbp, by=c("id", "visdate")) %>% 
  left_join(lung_func, by=c("id", "visdate")) %>% 
  left_join(renal_func, by=c("id", "visdate")) %>% 
  left_join(liver_fib, by=c("id", "visdate")) %>% 
  left_join(cesd, by=c("id", "visdate"))


# Define multimorbidity ---------------------------------------------------

# Someone flagged as having condition if 2 measurements within 1 year
dat4 <- dat3 %>%
  group_by(id) %>% 
  mutate(across(c(diab_idx, hbp_idx, lung_idx, renal_idx, liver_idx, dep_idx), 
                ~ ifelse(is.na(.x), 0, .x)),
         diabetes = as.numeric((cumsum(diab_idx + diabtx6m))>=1),
         hypertension = as.numeric((cumsum(hbp_idx + hbptx6m))>=1),
         lung_disease = as.numeric(cumsum(lung_idx)>=1),
         renal_disease = as.numeric(cumsum(renal_idx)>=1),
         depression = as.numeric((cumsum(dep_idx + deptx6m))>=1),
         liver_disease = as.numeric(cumsum(liver_idx)>=1)) %>% 
  select(-c(n_vis, yr, dx_date, poor_lung_func, renal_dysfunc, fibrosis, cirrhosis,
            hi_bp, cesd23))

# Due to stability of measurement requirement, drop everyone's first record
dat5 <- filter(dat4, duplicated(id))


# Build age 50 data -------------------------------------------------------

# Grab index visit
age50 <- filter(dat5, !duplicated(id, fromLast=T))

write_csv(age50, "../data/alive_cluster_setup_age50.csv")


# Build annual data -------------------------------------------------------

# Determine year participant had diagnosis
year.diagnosis <- dat5 %>% 
  group_by(id) %>% 
  mutate(yr = year(visdate),
         cum_hyper = cumsum(hypertension),
         year_hyper = ifelse(cum_hyper==1, yr, NA),
         cum_diab = cumsum(diabetes),
         year_diab = ifelse(cum_diab==1, yr, NA),
         cum_lung = cumsum(lung_disease),
         year_lung = ifelse(cum_lung==1, yr, NA),
         cum_renal = cumsum(renal_disease),
         year_renal = ifelse(cum_renal==1, yr, NA),
         cum_liver = cumsum(liver_disease),
         year_liver = ifelse(cum_liver==1, yr, NA),
         cum_dep = cumsum(depression),
         year_dep = ifelse(cum_dep==1, yr, NA),
         cum_hiv = cumsum(hiv),
         year_hiv = ifelse(cum_hiv==1, yr, NA),
         cum_cancer = cumsum(cancer),         
         year_cancer = ifelse(cum_cancer==1, yr, NA)) %>% 
  select(id, year_hyper, year_diab, year_lung, year_renal, 
         year_liver, year_dep, year_hiv, year_cancer) %>% 
  na.locf(na.rm=F) %>% 
  filter(!duplicated(id, fromLast=T))  
  
# Create one record per person-year
annual <- dat5 %>%
  ungroup() %>% 
  expand(id, yr=year(visdate)) %>% 
  left_join(year.diagnosis, by="id") %>% 
  mutate(hypertension = ifelse(is.na(year_hyper) | year_hyper>yr, 0, 1),
         diabetes = ifelse(is.na(year_diab) | year_diab>yr, 0, 1),
         lung_disease = ifelse(is.na(year_lung) | year_lung>yr, 0, 1),
         renal_disease = ifelse(is.na(year_renal) | year_renal>yr, 0, 1),
         liver_disease = ifelse(is.na(year_liver) | year_liver>yr, 0, 1),
         depression = ifelse(is.na(year_dep) | year_dep>yr, 0, 1),
         hiv = ifelse(is.na(year_hiv) | year_hiv>yr, 0, 1),
         cancer = ifelse(is.na(year_cancer) | year_cancer>yr, 0, 1)) %>% 
  select(id, yr, hypertension, diabetes, lung_disease, renal_disease, liver_disease,
         depression, hiv, cancer)

# Remove data after someone was censored
last.visit <- dat5 %>% 
  filter(!duplicated(id, fromLast=T)) %>% 
  select(id, visdate)

annual2 <- left_join(annual, last.visit, by="id") %>% 
  filter(yr <= year(visdate))

# Remove data before someone entered data



# Output data -------------------------------------------------------------

write_csv(dat6, "../data/alive_cluster_setup_annual.csv")

