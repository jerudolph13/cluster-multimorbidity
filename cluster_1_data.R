
###########################################################################
#
# Project: Multimorbidity clusters
#
# Purpose: Set up data
#
# Author: Jacqueline Rudolph
#
# Last Update: 06 Feb 2023
#
###########################################################################

library("tidyverse")
library("lubridate")


# Read in data ------------------------------------------------------------

base <- read_csv("../data/alive_base.csv") %>% 
  select(id, beduchs, enroll_period)

socdem <- read_csv("../data/alive_socdem.csv") %>% 
  filter(year(visdate) %in% 2007:2019) %>% 
  select(id, visdate, black, age, m0f1, anyhomels, inclt5k, work)

comorbid <- read_csv("../data/alive_comorbid.csv") %>% 
  filter(year(visdate) %in% 2007:2019) %>% 
  select(id, visdate, hbptx6m, bpsys, bpdias, diabtx6m, hba1c,
         fev_fvc, gfr, uprt_crt, fbscan)

mental <- read_csv("../data/alive_mental.csv") %>% 
  filter(year(visdate) %in% 2007:2019) %>% 
  select(id, visdate, deptx6m, cesdtot)

drug <- read_csv("../data/alive_drug.csv") %>% 
  filter(year(visdate) %in% 2007:2019) %>% 
  select(id, visdate, alcheavy, alcmoderate, curuser, anntivwomj, cigyn)

hiv <- read_csv("../data/alive_hiv.csv") %>% 
  filter(year(visdate) %in% 2007:2019) %>% 
  select(id, visdate, hiv, hcvvis)

care <- read_csv(file="../data/alive_care.csv") %>% 
  filter(year(visdate) %in% 2007:2019) %>% 
  select(id, visdate, primcare, insur6m)


# Manage data -------------------------------------------------------------

# Grab baseline demographic variables
socdem.base <- socdem %>% 
  filter(!duplicated(id)) %>% 
  select(id, m0f1, black)

# Merge data 
comorbid2 <- comorbid %>% 
  left_join(base, by="id") %>% 
  left_join(socdem.base, by="id") %>% 
  left_join(mental, by=c("id", "visdate")) %>% 
  left_join(hiv, by=c("id", "visdate")) %>% 
  left_join(select(socdem, id, age, visdate, work, inclt5k, anyhomels), by=c("id", "visdate")) %>% 
  left_join(care, by=c("id", "visdate")) %>% 
  left_join(drug, by=c("id", "visdate")) %>% 
  group_by(id) %>% 
  mutate(n = 1,
         n_vis = cumsum(n)) %>% 
  select(-n)
  
# Find index visit
index <- comorbid2 %>% 
  group_by(id) %>% 
  mutate(first = as.numeric(!duplicated(id))) %>% 
  filter(floor(age)>=50) %>% 
  select(id, visdate, age, first, n_vis)

  # Definition 1: first visit aged 50, if not first visit observed
  def1 <- filter(index, first!=1 & floor(age)==50) %>% 
    filter(!duplicated(id)) %>% 
    select(id, visdate, n_vis)
  
  # Definition 2: if person aged 50 at first visit observed, use 2nd visit within 365 days
  idx <- index$id[index$first==1 & floor(index$age)==50]
  def2 <- filter(index, id %in% idx) %>%
    filter(!(id %in% def1$id)) %>% 
    mutate(t_between = (lag(visdate) %--% visdate)/days(1)) %>% 
    filter(!is.na(t_between)) %>% 
    filter(!duplicated(id) & t_between<365) %>% 
    select(id, visdate, n_vis)
  
index2 <- bind_rows(def1, def2) %>% 
  mutate(flag = 1) %>% 
  rename(index_vis = n_vis)

# Merge back and remove visits after index date
comorbid3 <- comorbid2 %>% 
  filter(id %in% index2$id) %>% 
  left_join(select(index2, -index_vis), by=c("id", "visdate")) %>% 
  left_join(select(index2, index_vis), by=c("id")) %>% 
  group_by(id) %>% 
  mutate(flag = ifelse(is.na(flag), 0, flag),
         cum_flag = cumsum(cumsum(flag))) %>% 
  filter(cum_flag<=1)

# Inclusion criteria: 
  # Must have 1 visit within 365 days of index date
include <- comorbid3 %>% 
  group_by(id) %>% 
  mutate(t_between = (lag(visdate) %--% visdate)/days(1),
         rule = as.numeric(flag==1 & t_between<365))
ids <- filter(include, rule==1)

comorbid4 <- filter(comorbid3, id %in% ids$id) %>% 
  mutate(delta_vis = index_vis - n_vis) %>% 
  select(-c(flag, cum_flag, index_vis, delta_vis))


# Output data -------------------------------------------------------------

write_csv(comorbid4, "../data/alive_cluster_long.csv")


  