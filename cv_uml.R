
################################################################################
#
# Project: Multimorbidity clusters
#
# Purpose: v-fold cross-validation for unsupervised machine learning
#
# Last Update: 13 Feb 2022
#
################################################################################


# Cross-validation --------------------------------------------------------

uml_cv <- function(data, algorithm, k, nsplits) {
  set.seed(k)
  
  if (nsplits==1) {
    
    data.num <- data %>% 
      mutate(across(everything(), ~ as.numeric(.x)))
    
    # Run algorithm in full data
    if (algorithm=="kmeans") {
      data.res <- kmeans(data, centers=k, nstart=50)
      data.cluster <- data.res$cluster
    } else if (algorithm=="hierarchical") {
      diss.matrix <- daisy(data, metric="gower")
      data.res <- hclust(diss.matrix, method = "complete")
      data.cluster <- cutree(data.res, k=k)
    } else if (algorithm=="pam") {
      diss.matrix <- daisy(data, metric="gower")
      data.res <- pam(diss.matrix, k=k, diss=T)
      data.cluster <- data.res$clustering
    } else if (algorithm=="pdq") {
      # This could be generalized to figure out the variable types
      data.res <- PDQ(x=data, k=k, ini="random", dist="gower", bin=1:ncol(data))
      data.cluster <- data.res$label
    }
    
    res.k <- data.frame(k = k,
                        ps = "Undefined",
                        ch = NA,
                        si = NA)
    
    res.k$ch <- calinhara(data.num, data.cluster) # LIMITATION: CH Index seems to require numeric data
    if (k==1) {
      res.k$si <- NA
    } else {
      res.k$si <- mean(silhouette(data.cluster, daisy(data, metric="gower"))[ , 3])
    }
    
  } else if (nsplits>1) {
    
    # Split the data
    split <- sample(rep(1:nsplits, ceiling(nrow(data)/nsplits))[1:nrow(data)])
    
    # Data frame to hold results across v
    res.v <- data.frame(v=1:nsplits, 
                        ps=rep(NA, nsplits), 
                        ch=rep(NA, nsplits), 
                        si=rep(NA, nsplits))
    
    for (v in 1:nsplits) {
      
      val <- data[split==v, ]
      val.num <- val %>% 
        mutate(across(everything(), ~ as.numeric(.x))) # Make numeric for computation of means
      
      train <- data[split!=v, ]
      train.num <- train %>% 
        mutate(across(everything(), ~ as.numeric(.x))) # Make numeric for computation of means
      
      # Run algorithm in validation data
      if (algorithm=="kmeans") {
        val.res <- kmeans(val, centers=k, nstart=50)
        val.cluster <- val.res$cluster
      } else if (algorithm=="hierarchical") {
        diss.matrix <- daisy(val, metric="gower")
        val.res <- hclust(diss.matrix, method="complete")
        val.cluster <- cutree(val.res, k=k)
      } else if (algorithm=="pam") {
        diss.matrix <- daisy(val, metric="gower")
        val.res <- pam(diss.matrix, k=k, diss=T)
        val.cluster <- val.res$clustering
      } else if (algorithm=="pdq") {
        # This could be generalized to figure out the variable types
        val.res <- PDQ(x=val, k=k, ini="random", dist="gower", bin=1:ncol(val))
        val.cluster <- val.res$label
      }
      
      # Run algorithm in training data
      if (algorithm=="kmeans") {
        train.res <- kmeans(train, centers=k, nstart=50)
        train.cluster <- train.res$cluster
      } else if (algorithm=="hierarchical") {
        diss.matrix <- daisy(train, metric="gower")
        train.res <- hclust(diss.matrix, method = "complete")
        train.cluster <- cutree(train.res, k=k)
      } else if (algorithm=="pam") {
        diss.matrix <- daisy(train, metric="gower")
        train.res <- pam(diss.matrix, k=k, diss=T)
        train.cluster <- train.res$clustering
      } else if (algorithm=="pdq") {
        # This could be generalized to figure out the variable types
        train.res <- PDQ(x=train, k=k, ini="random", dist="gower", bin=1:ncol(train))
        train.cluster <- train.res$label
      }
      
      # Predict what training data group the validation ids would be assigned to, 
      #    based on minimizing distance to group centroid
      if (algorithm=="kmeans") {
        
        pred <- cl_predict(train.res, val)
        
      } else if (algorithm %in% c("hierarchical", "pam", "pdq")) {
        
        # Where are the training data centroids?
        train.centroid <- matrix(nrow=k, ncol=ncol(train))
        for (r in 1:k) { # For each group
          for (s in 1:ncol(train)) { # For each variable
            train.centroid[r, s] <- mean(train.num[train.cluster==r, ][[s]])
          }
        }
        
        # What is the sum of squared errors?
        sse <- matrix(nrow=nrow(val), ncol=k)
        for (t in 1:nrow(val)) {
          for (u in 1:k) {
            sse[t, u] <- sum((as.matrix(val.num[t, ]) - train.centroid[u, ])^2)
          }
        }
        pred <- apply(sse, 1, which.min) 
        
      } 
      
      # Compare assignments for observations in validation data
      val2 <- bind_cols(val, val_res=val.cluster, pred_res=pred)
      
      ps <- data.frame(k=1:k, res = rep(NA, k))
      for (i in 1:k) {
        # For each group in the validation subset
        subset <- filter(val2, val_res==i)
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
      
      res.v$ps[v] <- min(ps$res, na.rm=T) 
      res.v$ch[v] <- calinhara(val.num, pred) # LIMITATION: CH Index seems to require numeric data
      if (k==1) {
        res.v$si[v] <- NA
      } else {
        res.v$si[v] <- mean(silhouette(pred, daisy(val, metric="gower"))[ , 3])
      }
      
    }
    
    res.k <- res.v %>% 
      summarize(mean_ps = mean(ps),
                mean_ch = mean(ch),
                min_ch = min(ch),
                mean_si = mean(si),
                min_si = min(si)) %>% 
      mutate(k = k)
  }
 
return(res.k)

}

