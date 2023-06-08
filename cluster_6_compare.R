
###########################################################################
#
# Project: Multimorbidity clusters in ALIVE
#
# Purpose: Compare clustering approaches
#
# Author: Jacqueline Rudolph
#
# Last Update: 06 June 2023
#
###########################################################################


packages <- c("tidyverse", "lubridate", "caret", "ranger", "glmnet", "randomForest", "randomForestSRC",
              "RhpcBLASctl", "e1071", "xgboost", "SuperLearner", "cluster", "fpc")
for (package in packages) {
  library(package, character.only=T)
}

# Which model?
model <- "covariates alone" # {cluster alone, with covariates, covariates alone}


# What clusters do we want to describe?
algorithm <- "ensemble" # {"pam", "hier", "pdq", "ensemble"}
k <- 4


# Read in data ------------------------------------------------------------

res <- read_csv(file="../results/res_cluster_all.csv")

qual <- read_csv(file="../data/alive_qual.csv") %>% 
  select(id, visdate, mosmhs) # or mosphs

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

dat <- res %>% 
  left_join(qual, by=c("id", "visdate")) %>% 
  mutate(delta = mosmhs, # or mosphs
         delta_bin = as.numeric(mosmhs<50))

if (model=="cluster alone") {
  dat2 <- dat %>%
    mutate(cluster = factor(cluster)) %>% 
    select(cluster, delta) %>%
    na.omit()
} else if (model=="with covariates") {
  dat2 <- dat %>%
    mutate(cluster = factor(cluster), 
           across(c(black, m0f1, enroll_period, inclt5k, beduchs, curuser, alcheavy, cigyn), ~ factor(.x))) %>% 
    select(cluster, delta, black, m0f1, enroll_period, inclt5k, beduchs, curuser, alcheavy, cigyn) %>%
    na.omit()
} else if (model=="covariates alone") {
  dat2 <- dat %>%
    mutate(across(c(black, m0f1, enroll_period, inclt5k, beduchs, curuser, alcheavy, cigyn), ~ factor(.x))) %>% 
    select(delta, black, m0f1, enroll_period, inclt5k, beduchs, curuser, alcheavy, cigyn) %>%
    na.omit()
}



# Split the data ----------------------------------------------------------

set.seed(123)

# Split the data
training.indices <- createDataPartition(dat2$delta, p=0.7, list=F)
train.data <- dat2[training.indices, ]
test.data <- dat2[-training.indices, ]


# Prediction model --------------------------------------------------------

# SuperLearner
Y <- train.data$delta
X <- select(train.data, -delta)
learn.xgboost <- create.Learner("SL.xgboost", tune = list(ntrees = c(10, 20),
                                                     max_depth = 1:2,
                                                     shrinkage = c(0.001, 0.01)), 
                           detailed_names = TRUE, name_prefix = "xgb")
learn.ranger <- create.Learner("SL.ranger", tune = list(num.trees = c(100, 500, 1000), 
                                                        mtry = seq(1, ncol(dat2)-1, 1)),
                               detailed_names = TRUE, name_prefix = "ranger")
learn.glmnet <- create.Learner("SL.glmnet", tune = list(alpha = c(0.0, 0.001, 0.0025, 0.005, 0.0075, 0.01, 
                                                                  0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 
                                                                  0.75, 1.0)),
                               detailed_names = TRUE, name_prefix = "glmnet")

library <- c("SL.mean", "SL.glm", learn.glmnet$names, learn.ranger$names, learn.xgboost$names)
sl <- SuperLearner(Y=Y, X=X, family=gaussian(), SL.library=library, method = "method.NNLS")
sl

# MSE in the training data
train.pred <- sl$SL.predict
mean((train.pred - train.data$delta)^2)

# MSE in the testing data
test.pred <- predict(sl, select(test.data, -delta), onlySL = TRUE)
mean((test.pred$pred - test.data$delta)^2)
  
  
# AUC if outcome is binary
# pred_rocr <- ROCR::prediction(pred$pred, test.data$delta)
# auc <- ROCR::performance(pred_rocr, measure = "auc", x.measure = "cutoff")@y.values[[1]]
# auc
# pROC::ci.auc(as.numeric(test.data$delta), as.numeric(pred$pred))


