
###########################################################################
#
# Project: Multimorbidity clusters in ALIVE
#
# Purpose: Make figures
#
# Author: Jacqueline Rudolph
#
# Last Update: 16 Feb 2023
#
###########################################################################

packages <- c("tidyverse", "gmodels", "gifski", "gganimate", 
              "magick", "ggpubr", "transformr", "patchwork")
for (package in packages) {
  library(package, character.only=T)
}

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


# Manage data -------------------------------------------------------------

# Frequency of clusters
CrossTable(x=res$cluster)

# Median number of conditions by cluster
medians <- res %>% 
  group_by(cluster) %>% 
  mutate(n_cond = hypertension + diabetes + lung_disease + renal_disease +
         liver_disease + depression + cancer + hiv) %>% 
  summarize(med_cond = median(n_cond)) %>% 
  mutate(cluster = factor(cluster))

# Disease distribution overall
ranks1 <- res %>% 
  summarize(hyper = mean(hypertension),
            depression = mean(depression),
            hiv = mean(hiv),
            renal = mean(renal_disease),  
            diab = mean(diabetes),
            lung = mean(lung_disease),
            liver = mean(liver_disease),
            cancer = mean(cancer)) %>% 
  pivot_longer(everything(), names_to="condition", values_to="overall_prevalence")

# Disease distribution by cluster
ranks2 <- res %>% 
  group_by(cluster) %>% 
  summarize(hyper = mean(hypertension),
            depression = mean(depression),
            hiv = mean(hiv),
            renal = mean(renal_disease),  
            diab = mean(diabetes),
            lung = mean(lung_disease),
            liver = mean(liver_disease),
            cancer = mean(cancer)) %>% 
  pivot_longer(!cluster, names_to="condition", values_to="prevalence")

diseases <- c("Cancer", "Diabetes", "Liver Disease", "Lung Disease",   
              "Renal Disease", "HIV", "Depression", "Hypertension")
ranks <- ranks2 %>% 
  left_join(ranks1, by="condition") %>% 
  mutate(cluster = factor(cluster),
         diff_mean = prevalence - overall_prevalence,
         condition = factor(condition, 
                            levels=c("cancer", "diab", "liver", "lung", "renal", "hiv", 
                                     "depression", "hyper"),
                            labels=diseases),
         condition_no = recode(condition, 
                               "Cancer" = 1,
                               "Diabetes" = 2, 
                               "Liver Disease" = 3,
                               "Lung Disease" = 4,
                               "Renal Disease" = 5,
                               "HIV" = 6,
                               "Depression" = 7,
                               "Hypertension" = 8)) %>% 
  left_join(medians, by="cluster")


# Plot theme --------------------------------------------------------------

thm1 <- theme_classic() +
  theme(
    # Format plot title
    plot.title = element_text(family="Helvetica", size=16, color="black"),
    
    # Format axes
    axis.title = element_text(family="Helvetica", size=16, color="black"),
    axis.text.y = element_text(family="Helvetica", size=14, color="black"),
    axis.text.x = element_text(family="Helvetica", size=14, color="black"),
    axis.line = element_line(size=0.75),
    axis.ticks = element_line(size=0.75),
    
    # Format legend
    legend.text = element_text(family="Helvetica", size=14, color="black", margin=margin(t=0.25,b=0.25, unit="lines")),
    legend.title = element_text(family="Helvetica", size=16, color="black"),
    legend.title.align = 0.5,
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.direction = "horizontal",
    
    # Format tags
    plot.tag = element_text(family="Helvetica", size=16, color="black"),
    
    # Add space around plot
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    
    #Format facet title
    strip.text = element_text(family="Helvetica", size=14, color="black"),
    strip.background = element_blank(),
    # panel.border = element_rect(size=0.75, fill=NA),
    panel.spacing = unit(30, "pt")
  )

# Disease distribution by cluster -----------------------------------------

fig <- rename(ranks, Cluster=cluster)

# jpeg(paste0("../figures/diseases_", algorithm, k,".jpeg"), 
#      height=5*(as.numeric(k>3) + 1), 
#      width=ifelse(k>3, 12, 4*k), 
#      units="in", res=300)
jpeg(paste0("../figures/diseases_", algorithm, k,".jpeg"), 
     height=5, 
     width=4*k, 
     units="in", res=300)
ggplot(data=fig, aes(y=condition, x=prevalence*100, fill=condition)) + thm1 +
  labs(x="\nPrevalence at age 50", y="", ) +
  geom_bar(stat="identity", show.legend=F) +
  scale_x_continuous(expand=c(0, 0), limits=c(0, 100)) +
  # scale_fill_grey() +
  facet_wrap(~Cluster, ncol=5, nrow=2, labeller=label_both)
dev.off()


