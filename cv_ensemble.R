
################################################################################
#
# Project: Multimorbidity clusters
#
# Purpose: v-fold cross-validation for ensemble unsupervised machine learning
#
# Last Update: 06 Jun 2022
#
################################################################################

  

# Cross-validation --------------------------------------------------------

uml_cv <- function(obs, clust, alpha1, alpha2, nsplits) {
  
  # Split the data
  set.seed(123)
  split <- sample(rep(1:nsplits, ceiling(nrow(obs)/nsplits))[1:nrow(obs)])
  
  # Data frame to hold results across v
  res.v <- data.frame(v=1:nsplits, 
                      ps=rep(NA, nsplits), 
                      ch=rep(NA, nsplits), 
                      si=rep(NA, nsplits))
  
  for (v in 1:nsplits) {
    
    val.obs <- obs[split==v, ]
      val.num <- val.obs %>% 
        mutate(across(everything(), ~ as.numeric(.x)))
    val.clust <- clust[split==v, ]
    
    train.obs <- obs[split!=v, ]
      train.num <- train.obs %>% 
        mutate(across(everything(), ~ as.numeric(.x)))
    train.clust <- clust[split!=v, ]

    # Run algorithm in validation data
    res1 <- cl.sim(tmp=val.clust, bin.ind=T, wt=NULL)
    res2 <- merge.clust(cls.res=res1, alpha=alpha1)
      # Calculate membership similarity
      member <- memb.sim(zdat=res2)
      # Assign cluster based on membership similarity
      new.clust <- assign.clust(zdat=member, alpha2=alpha2, wt=rep(1, nrow(member)))
      val.cluster <- new.clust[[1]]
    
    # Run algorithm in training data
    res1 <- cl.sim(tmp=train.clust, bin.ind=T, wt=NULL)
    res2 <- merge.clust(cls.res=res1, alpha=alpha1)
      # Calculate membership similarity
      member <- memb.sim(zdat=res2)
      # Assign cluster based on membership similarity
      new.clust <- assign.clust(zdat=member, alpha2=alpha2, wt=rep(1, nrow(member)))
      train.cluster <- new.clust[[1]]
      
    # Predict what training data group the validation ids would be assigned to, 
    #    based on minimizing distance to group centroid
    # Where are the training data centroids?
    k <- as.numeric(length(table(train.cluster)))
    train.centroid <- matrix(nrow=k, ncol=ncol(train.obs))
    for (r in 1:k) { # For each group
      for (s in 1:ncol(train.obs)) { # For each variable
        train.centroid[r, s] <- mean(train.num[train.cluster$cluster==r, ][[s]])
      }
    }
      
    # What is the sum of squared errors?
    sse <- matrix(nrow=nrow(val.obs), ncol=k)
    for (t in 1:nrow(val.obs)) {
      for (u in 1:k) {
        sse[t, u] <- sum((as.matrix(val.num[t, ]) - train.centroid[u, ])^2)
      }
    }
    pred <- apply(sse, 1, which.min) 
    
    # Compare assignments for observations in validation data
    val2 <- bind_cols(val.obs, val_res=val.cluster[[1]], pred_res=pred)
    
    k2 <- as.numeric(length(table(val.cluster))) # Just in case # clusters differs
    ps <- data.frame(k=1:k2, res = rep(NA, k2))
    for (i in 1:k2) {
      # For each group in the validation subset
      subset <- filter(val2, val_res==i)
      
      if (nrow(subset)==0) {
        ps$res[i] <- NA
      } else {
        D <- matrix(nrow=nrow(subset), ncol=nrow(subset))
        
        # See if individuals are in same group in training data
        for (i1 in 1:nrow(subset)) {
          for (i2 in 1:nrow(subset)) {
            D[i1, i2] <- as.numeric(subset$pred_res[i1]==subset$pred_res[i2])
          }
        }
        
        matches <- sum(D[row(D)!=col(D)])
        ps$res[i] <- (1/(nrow(subset)*(nrow(subset)-1)))*matches
      }
    }
    
    res.v$ps[v] <- min(ps$res, na.rm=T)
    res.v$ch[v] <- calinhara(val.num, pred)
    # For SI, handle case where all obs assinged to same cluster
    si <- silhouette(pred, daisy(val.obs, metric="gower"))
    res.v$si[v] <- ifelse(is.na(si), NA, mean(silhouette(pred, daisy(val.obs, metric="gower"))[ , 3]))

  }
  
  res.k <- res.v %>% 
    dplyr::summarize(mean_ps = mean(ps, na.rm=T),
              mean_ch = mean(ch, na.rm=T),
              min_ch = min(ch, na.rm=T),
              mean_si = mean(si, na.rm=T),
              min_si = min(si, na.rm=T)) %>% 
    mutate(alpha1 = alpha1,
           alpha2 = alpha2)
  
  return(res.k)
  
}

