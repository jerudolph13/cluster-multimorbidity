
###########################################################################
#
# Project: Multimorbidity clusters in ALIVE
#
# Purpose: Make tables
#
# Author: Jacqueline Rudolph
#
# Last Update: 06 June 2023
#
###########################################################################


library("tidyverse")
library("gmodels")

# What clusters do we want to describe?
algorithm <- "ensemble" # {"pam", "hier", "pdq", "ensemble"}
k <- 4


# Read in results ---------------------------------------------------------

res <- read_csv(file="../results/res_cluster_all.csv")

# Which clusters do we want to describe?
if (algorithm=="pam" & k==2) {
  res$cluster <- res$pam2
} else if (algorithm=="pam" & k==3) { 
  res$cluster <- res$pam3
} else if (algorithm=="pam" & k==5) { 
  res$cluster <- res$pam5
} else if (algorithm=="hier" & k==2) {
  res$cluster <- res$hier2
} else if (algorithm=="hier" & k==3) {
  res$cluster <- res$hier3
} else if (algorithm=="pdq" & k==2) {
  res$cluster <- res$pdq2
} else if (algorithm=="pdq" & k==5) {
  res$cluster <- res$pdq5
} else if (algorithm=="ensemble") {
  res$cluster <- res$ensemble
}


# Trends ------------------------------------------------------------------

# Cluster frequency
CrossTable(res$cluster, prop.r=T, prop.c=T, prop.t=T)
  

# Health conditions -------------------------------------------------------

# Distribution overall
cond.overall <- res %>%
  mutate(group = "Overall", 
         n_cond = hypertension + diabetes + lung_disease + renal_disease +
                  liver_disease + depression + cancer + hiv) %>% 
  group_by(group) %>% 
  summarize(n_cond = paste0(round(median(n_cond), 1), " (", 
                            round(quantile(n_cond, probs=0.25), 1), ", ", 
                            round(quantile(n_cond, probs=0.75), 1), ")"),
            hyper = paste0(sum(hypertension==1, na.rm=T), " (", round(100*sum(hypertension==1, na.rm=T)/n(), 1), ")"),
            depression = paste0(sum(depression==1, na.rm=T), " (", round(100*sum(depression==1, na.rm=T)/n(), 1), ")"),
            hiv = paste0(sum(hiv==1, na.rm=T), " (", round(100*sum(hiv==1, na.rm=T)/n(), 1), ")"),
            renal_disease = paste0(sum(renal_disease==1, na.rm=T), " (", round(100*sum(renal_disease==1, na.rm=T)/n(), 1), ")"),
            lung_disease = paste0(sum(lung_disease==1, na.rm=T), " (", round(100*sum(lung_disease==1, na.rm=T)/n(), 1), ")"),
            liver_disease = paste0(sum(liver_disease==1, na.rm=T), " (", round(100*sum(liver_disease==1, na.rm=T)/n(), 1), ")"),
            diab = paste0(sum(diabetes==1, na.rm=T), " (", round(100*sum(diabetes==1, na.rm=T)/n(), 1), ")"),
            cancer = paste0(sum(cancer==1, na.rm=T), " (", round(100*sum(cancer==1, na.rm=T)/n(), 1), ")"))

# Distribution by cluster
cond.cluster <- res %>%
  mutate(group = paste0("Cluster ", cluster),
         n_cond = hypertension + diabetes + lung_disease + renal_disease +
           liver_disease + depression + cancer + hiv) %>% 
  group_by(group) %>% 
  summarize(n_cond = paste0(round(median(n_cond), 1), " (", 
                            round(quantile(n_cond, probs=0.25), 1), ", ", 
                            round(quantile(n_cond, probs=0.75), 1), ")"),
            hyper = paste0(sum(hypertension==1, na.rm=T), " (", round(100*sum(hypertension==1, na.rm=T)/n(), 1), ")"),
            depression = paste0(sum(depression==1, na.rm=T), " (", round(100*sum(depression==1, na.rm=T)/n(), 1), ")"),
            hiv = paste0(sum(hiv==1, na.rm=T), " (", round(100*sum(hiv==1, na.rm=T)/n(), 1), ")"),
            renal_disease = paste0(sum(renal_disease==1, na.rm=T), " (", round(100*sum(renal_disease==1, na.rm=T)/n(), 1), ")"),
            lung_disease = paste0(sum(lung_disease==1, na.rm=T), " (", round(100*sum(lung_disease==1, na.rm=T)/n(), 1), ")"),
            liver_disease = paste0(sum(liver_disease==1, na.rm=T), " (", round(100*sum(liver_disease==1, na.rm=T)/n(), 1), ")"),
            diab = paste0(sum(diabetes==1, na.rm=T), " (", round(100*sum(diabetes==1, na.rm=T)/n(), 1), ")"),
            cancer = paste0(sum(cancer==1, na.rm=T), " (", round(100*sum(cancer==1, na.rm=T)/n(), 1), ")"))

# Transform into table
cond <- bind_rows(cond.overall, cond.cluster) %>% 
  pivot_longer(!group, names_to="Variable", values_to="val") %>% 
  pivot_wider(id_cols="Variable", names_from="group", values_from="val")
cond$Variable <- recode(cond$Variable, 
                        "n_cond" = "No. conditions, median (IQR)",
                        "hyper" = "Hypertension",
                        "depression" = "Depression",
                        "hiv" = "HIV",
                        "renal_disease" = "Renal disease",
                        "lung_disease" = "Lung disease",
                        "liver_disease" = "Liver disease",
                        "diab" = "Diabetes",
                        "cancer" = "Cancer")

# Output table
write_csv(cond, paste0("../results/", algorithm, k, "_conditions.csv"))


# Participant characteristics ---------------------------------------------

# Distribution overall
CrossTable(res$enroll_period, res$cluster, prop.r=F, prop.c=T, prop.t=F, prop.chisq=F)

char.overall <- res %>%
  mutate(group = "Overall") %>% 
  group_by(group) %>% 
  summarize(female = paste0(sum(m0f1==1, na.rm=T), " (", round(100*sum(m0f1==1, na.rm=T)/n(), 1), ")"),
            black = paste0(sum(black==1, na.rm=T), " (", round(100*sum(black==1, na.rm=T)/n(), 1), ")"),
            income = paste0(sum(inclt5k==1, na.rm=T), " (", round(100*sum(inclt5k==1, na.rm=T)/n(), 1), ")"),
            education = paste0(sum(beduchs==1, na.rm=T), " (", round(100*sum(beduchs==1, na.rm=T)/n(), 1), ")"),
            idu = paste0(sum(curuser==1, na.rm=T), " (", round(100*sum(curuser==1, na.rm=T)/n(), 1), ")"),
            alc = paste0(sum(alcheavy==1, na.rm=T), " (", round(100*sum(alcheavy==1, na.rm=T)/n(), 1), ")"),
            cig = paste0(sum(cigyn==1, na.rm=T), " (", round(100*sum(cigyn==1, na.rm=T)/n(), 1), ")"))

# Distribution by cluster
char.cluster <- res %>%
  mutate(group = paste0("Cluster ", cluster)) %>% 
  group_by(group) %>% 
  summarize(female = paste0(sum(m0f1==1, na.rm=T), " (", round(100*sum(m0f1==1, na.rm=T)/n(), 1), ")"),
            black = paste0(sum(black==1, na.rm=T), " (", round(100*sum(black==1, na.rm=T)/n(), 1), ")"),
            income = paste0(sum(inclt5k==1, na.rm=T), " (", round(100*sum(inclt5k==1, na.rm=T)/n(), 1), ")"),
            education = paste0(sum(beduchs==1, na.rm=T), " (", round(100*sum(beduchs==1, na.rm=T)/n(), 1), ")"),
            idu = paste0(sum(curuser==1, na.rm=T), " (", round(100*sum(curuser==1, na.rm=T)/n(), 1), ")"),
            alc = paste0(sum(alcheavy==1, na.rm=T), " (", round(100*sum(alcheavy==1, na.rm=T)/n(), 1), ")"),
            cig = paste0(sum(cigyn==1, na.rm=T), " (", round(100*sum(cigyn==1, na.rm=T)/n(), 1), ")"))

# Transform into table
char <- bind_rows(char.overall, char.cluster) %>% 
  pivot_longer(!group, names_to="Variable", values_to="val") %>% 
  pivot_wider(id_cols="Variable", names_from="group", values_from="val")
char$Variable <- recode(char$Variable, 
                        "female" = "Female",
                        "black" = "Black",
                        "income" = "Income <$5000", 
                        "education" = "High school education",
                        "idu" = "Any injection drug use",
                        "alc" = "Heavy alcohol use",
                        "cig" = "Cigarette use")

# Output table
write_csv(char, paste0("../results/", algorithm, k, "_characteristics.csv"))

