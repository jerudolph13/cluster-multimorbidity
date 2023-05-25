
###################################################################################################################
#
# Purpose: Ensemble clustering functions
#
# Author: Bryan Lau
#
# Last Update: 21 Mar 2022
#
###################################################################################################################


##############################################################################
# Adaptive Cluster Ensemble (ACE)
# relevant paper is:
#   Alqurashi and Wang (2019) International Journal of Machine Learning and Cybernetics 10:1227-1246
#   https://doi.org/10.1007/s13042-017-0756-7


# The ACE algorithm is as follows
#   1) transform all clustering algorithms into binary indicators of the clustering algorithm
#       and class membership within the algorithm
#   2) calculate cluster similarity
#   3) merge clusters based upon cluster similarity and set level of alpha
#   4) re-estimate cluster similarity
#   5) merge if still meets set level of alpha
#   6) if no more merging and you have Z clusters and want k number of clusters
#       A) if Z=k then stop here
#       B) if Z<k then alpha may be too low and can increase this threshold
#       c) if Z>k then use following options
#           i)   accept Z>k as given threshold there is no more similar clusters
#           ii)  lower alpha threshold as "there is no gold standard threshold"
#                 to allow merger until you get to k clusters
#           iii) cluster elimination
#                 - if each top k clusters have at least 1 certain object based upon
#                   current value of alpha2 then keep top k clusters and eliminate
#                   the rest. Update the cluster similarity and membership similarity
#                 - otherwise adapt alpha2 to be maximum membership similarity of the
#                   kth cluster 
#                 (NOTE: it isn't specified but ordering of the clusters
#                   could be done using "cluster certainty" which is just the mean
#                   of the membership similarity of the observations in cluster)
#   7) enforce hard clustering
#       A) identify any totally uncertain observations (note this may happen with
#           cluster elimination in step 6)
#           ##########################################################################
#           # my personal feeling is to accually do i) under step 6c but potentially
#           # 6cii would also acceptable. I'd rather not due elimination... but
#           #  perhaps this would be guided by cluster certainty
#           ##########################################################################
#       B) identify any totally certain or certain observations and assign these
#           observations to those clusters
#       c) measure quality of candidate cluster - authors suggest any measure
#           could be used "in principle". Authors use the variance which in this
#           case is the average square differences in the mean where the mean is
#           the average membership similarity of observations within a cluster

##############################################################################


##############################################################################
# Function to calculate similarity between two clusters across algorithms
#   or across clusters after doing merges in the ensemble
cluster.similarity <- function(cj,ck,wtval){
  # number of individuals that are in the intersection of the two clusters
  # need cj and ck for the clusters being compared
  #   note if weighted then these should no longer be binary indicators
  # wtval is being passed to function as n is no longer the number of 
  #   observations but the sum of the weights
  
  cj.ck  <- length(cj[cj==1 & ck==1])
  cjn <- sum(cj)
  ckn <- sum(ck)
  CM <- cj.ck/sqrt(cjn*ckn)
  n <- sum(wtval)
  
  # essentially function uses the norm and inner products
  #similarity <- (n*CM - sqrt(cjn*ckn))/(sqrt((n-cjn)*(n-ckn)))
  #similarity2 <- (cj.ck - (cjn*ckn)/n)/sqrt(cjn*ckn*(1-cjn/n)*(1-ckn/n))
  similarity3 <- n*(cj%*%ck)/sqrt((cj%*%cj) * (ck%*%ck)) - sqrt((cj%*%cj) * (ck%*%ck))
  similarity3 <- similarity3/sqrt((n-cj%*%cj)*(n-ck%*%ck))
  #return(c(similarity,similarity2,similarity3))
  return(c(similarity3))
}


##############################################################################
# function to transform from cluster membership to binary indicators

member.t <- function(dat){
  # dat is dataframe of cluster algorithms and class membership
  
  # internal function returns a list
  small.fun <- function(x){
    res <- matrix(NA,nrow=length(x),ncol=max(x))
    for(val in unique(x)){
      res[,val] <- ifelse(x==val,1,0)
    }
    res <- as.data.frame(res)
    return(res)
  }
  member.list <- apply(dat,2,small.fun)
  for(mval in c(1:length(member.list))){
    vv <- length(names(member.list[[mval]]))
    names(member.list[[mval]]) <- paste(names(member.list)[mval],c(1:vv),sep=".")
  }
  return(member.list)  
}


##############################################################################
# One main function
# Function takes:
#   1) dataframe of observations and class membership across algorithms
#      to determine similarity of clusters across algorithms
#   2) or results from merging based upon already running 1 to see if clusters
#      should be merged even more so

cl.sim <- function(tmp,bin.ind=T,wt=NULL){
  #wrapper around transformation and cluster.similarity
  # tmp is the dataframe of observations and class membership across algorithms
  
  # if function is not passed a weight then it is weight of 1
  if(is.null(wt) & bin.ind==T){
    wt <- rep(1, nrow(tmp))
  } else if (is.null(wt) & bin.ind==F) {
    wt <- rep(1, nrow(tmp[[1]]))
  }
  
  # transform to binary indicator list
  if(bin.ind==T) ind.list <- member.t(tmp)
  
  # if running on output of merge.clust but no merges occurred,
  # all lists are binary
  if(bin.ind==F & !any(unlist(tmp)>1)) {
    ind.list <- tmp
    
    # rename tmp and ind.list
    names(tmp) <- paste("c", c(1:length(tmp)), sep="")
    names(ind.list) <- paste("c", c(1:length(ind.list)), sep="")
    for(val in c(1:length(ind.list))) {
      names(tmp[[val]]) <- paste(names(tmp)[val], 1, sep=".")
      names(ind.list[[val]]) <- paste(names(ind.list)[val],1,sep=".")
    }
  # if running on output of merge.clust and merges occurred,
  #   create indicators of whether observation is within cluster or not
  } else if(bin.ind==F & any(unlist(tmp)>1)) {
    ind.list <- lapply(tmp, function(x) data.frame(as.numeric(x>0)))
    
    # rename tmp and ind.list
    names(tmp) <- paste("c", c(1:length(tmp)), sep="")
    names(ind.list) <- paste("c", c(1:length(ind.list)), sep="")
    for(val in c(1:length(ind.list))) {
      names(tmp[[val]]) <- paste(names(tmp)[val], 1, sep=".")
      names(ind.list[[val]]) <- paste(names(ind.list)[val],1,sep=".")
    }
  }
  
  # to allow for sample weights - make sure that tmp when bin.ind=F is already weighted
  ind.list <- lapply(ind.list, function(x) x*wt)
  
  # create the similarity matrix dimension
  col.n <- unlist(lapply(ind.list, ncol))
  sumcol <- cumsum(col.n)
  s.mat <- as.data.frame(matrix(NA, nrow=sum(col.n), ncol=sum(col.n)))
  listval <- c(1:length(ind.list))
  
  tmpfun <- function(xvec,ymat,xnm,ynm,wt.vec) {
    cvec <- rep(NA, ncol(ymat))
    for(qq in c(1:ncol(ymat))){
      cvec[qq] <- cluster.similarity(xvec ,ymat[,qq], wtval=wt.vec)
    }
    names(cvec) <- paste(xnm, ynm, c(1:ncol(ymat)), sep="")
    return(cvec)
  }
  
  rowstrt <- c(0, sumcol[-length(sumcol)]) + 1
  rowend <- sumcol
  colstrt <- sumcol[-length(sumcol)] + 1
  colend <- sumcol[-1]
  
  for(val in c(listval[-max(listval)])){
    for(val2 in c((val+1):max(listval))){
      tmpres <- apply(ind.list[[val]], 2, FUN=tmpfun, ymat=ind.list[[val2]],
                      xnm=names(ind.list[val]), ynm=names(ind.list[val2]), wt.vec=wt)
      lx <- c(rowstrt[val]:rowend[val])
      ly <- c(colstrt[val2-1]:colend[val2-1])
      s.mat[lx, ly] <- t(tmpres)
      dimnames(s.mat)[[2]][ly] <- names(ind.list[[val2]])
      dimnames(s.mat)[[1]][lx] <- names(ind.list[[val]])
    }
  }
  
  # because of how subset names - the first set of column indicators are not named
  # and the last subset of row indicators are not named
  dimnames(s.mat)[[1]][(ncol(s.mat)-col.n[length(col.n)]+1):ncol(s.mat)] <- 
    names(ind.list[[length(ind.list)]])
  dimnames(s.mat)[[2]][1:col.n[1]] <- names(ind.list[[1]])
  
  if(bin.ind==T) return(list(simularity=s.mat, binary.ind=ind.list))
  # removed condition any(unlist(tmp)>1)
  # when this condition is not met, cl.merge==binary.ind
  if(bin.ind==F) return(list(simularity=s.mat, cl.merge=tmp, binary.ind=ind.list))
}


##############################################################################
# Merge clusters based upon cluster similarity
merge.clust <- function(cls.res,alpha=0.8){
  #cls.res is a result from cl.sim function with the cluster similarity as element 1
  # and class indicator list as element 2
  
  sim.mat <- cls.res[[1]]
  
  # make the class indicator list into a dataframe
  ind.list <- bind_cols(cls.res[[2]])
  
  #make sim.mat into a symmetrical matrix
  sim.mat[lower.tri(sim.mat)] <- t(sim.mat)[lower.tri(sim.mat)]
  
  # which is above threshold
  gt.mat <- sim.mat>=alpha
  # # are any columns all FALSE?
  # f.res <- apply(gt.mat,2,table)
  # f.res <- unlist(lapply(f.res,length))
  # f.res <- names(f.res)[f.res==1]
  
  m.list <- vector(mode="list",length=ncol(sim.mat))
  # if column is all below alpha then will be it's own cluster for next step
  for(val in c(1:nrow(gt.mat))){
    curr.row <- rownames(gt.mat)[val]
    if(!any(unlist(lapply(m.list,function(x) curr.row %in% x)))){
      m.list[[val]] <- c(curr.row,colnames(gt.mat)[gt.mat[val,]==T & !is.na(gt.mat[val,])])
    }
  }
  m.list <- m.list[!sapply(m.list,is.null)]
  
  # combine clusters
  quickfun <- function(xvec,i.list){
    sum.these <- names(i.list) %in% xvec
    if(table(sum.these)[2]>=2){
      new.vec <- rowSums(i.list[,sum.these])
      new.vec <- as.data.frame(new.vec)
    } 
    if(table(sum.these)[2]<2) {
      tmp <- as.data.frame(i.list[,sum.these])
      names(tmp) <- xvec
      new.vec <- tmp
    }
    return(new.vec)
  }
  new.ind <- lapply(m.list,quickfun,i.list=ind.list)
  return(new.ind)
}


##############################################################################
# Membership similarity
# This function is seeing how similar an observation is to each of the clusters
# this is actually the sum of the indicator of the merged clusters divided by
# the a sum of the indicator overall clusters. This can actually be done on 
# the the resulting list from merge.clust

memb.sim <- function(zdat){
  # pass this function the resulting merged clusters from 
  # merge.clust function
  res <- bind_cols(zdat)
  res <- apply(res,1,function(x) x/sum(x))
  res <- data.frame(t(res))
  names(res) <- paste("c",c(1:ncol(res)),sep="")
  rownames(res) <- c(1:nrow(res)) 
  return(res)
}


##############################################################################
# enforcing Hard clustering

# assign observations to cluster above alpha2 threshold
#  let's say you have so many clusters that you might lower threshold below
#   0.5... I am not sure why one would do this, but in this case there is
#   the possibility that multiple clusters would be above threshold
#   going to set it to go with the cluster with the max value above the threshold

assign.clust <- function(zdat,alpha2=0.5,wts=NULL){
  # take the data frame of MEMBERSHIP SIMILARITY and pass to apply
  # which if any is above threshold of alpha2 is considered certain and assigned to
  # cluster
  
  # zdat - membership similarity from function memb.sim
  # alpha2 - threshold for determining certain and uncertain observations
  # wts - weights for sampling/IPW/etc.
  
  if(is.null(wts)) wts <- rep(1,nrow(zdat))
  
  res <- apply(zdat,1,function(x) ifelse(any(x>alpha2),c(1:ncol(zdat))[x>alpha2 & x==max(x)],NA))
  
  # if any uncertain observations (i.e. membership similarity below alpha2)
  # calculate cluster quality of those assigned
  #   cluster quality is essentially the variance of the cluster (see below)
  cqual <- clust.qual(zdat,res)
  
  #for uncertain observations, assign to each cluster that they were associated with
  # and calculate new cluster quality assign this observation to the cluster which
  # results in the least impact on quality
  #So starting with identifying the unassigned observations create a for loop that
  # cycles through these observations
  for(val in c(1:nrow(zdat))[is.na(res)]){
    # which clusters can they be assigned to?
    aclust <- zdat[val,]>0
    newqual <- rep(NA,ncol(zdat))
    if(!is.nan(zdat[val,1])){
      
      # For each of its possible clusters, save what the new cluster quality would be in newqual
      for(val2 in c(1:ncol(zdat))[aclust]){
        newres <- res
        newres[val] <- val2
        newqual[val2] <- clust.qual(zdat,newres,wtval=wts)[val2] 
      }
      
      # for the possible clusters, did they have a non-missing quality value?
      # if so, compare new quality to old quality
      if(!any(is.na(cqual[aclust]))){
        difqual <- newqual[aclust]-cqual[aclust]
        fclust <- c(1:ncol(zdat))[difqual==min(difqual,na.rm = T) & !is.na(difqual)]
        res[val] <- fclust
      }
      
      # so if first observation being added to cluster than initial cqual is NA
      if(!any(is.na(cqual[aclust]))){
        difqual <- newqual[aclust]-cqual[aclust]
        difqual[is.na(difqual)] <- Inf # setting at Inf because if first observation being added to 
        # cluster then all individuals were <alpha should add to 
        # the clusters that had some indication of membership belonging
        if(!all(difqual==Inf)){
          fclust <- c(1:ncol(zdat))[difqual==min(difqual,na.rm = T) & !is.na(difqual)]
        }
        if(all(difqual==Inf)) {
          # if all of the potential clusters have not been assigned yet
          #   then just do a random one
          fclust <- sample(dimnames(aclust)[[2]][aclust==T],1)
        }
      }
      
      # I'M NOT SURE IF BOTH OF THE ABOVE NEED TO BE THERE
      
      if(any(is.na(cqual[aclust]))){
        difqual <- newqual[aclust]-cqual[aclust]
        difqual[is.na(difqual)] <- Inf # setting at Inf because if first observation being added to 
        # cluster then all individuals were <alpha should add to 
        # the clusters that had some indication of membership belonging
        if(!all(difqual==Inf)){
          fclust <- c(1:ncol(zdat))[difqual==min(difqual,na.rm = T) & !is.na(difqual)]
        }
        if(all(difqual==Inf)) {
          # if all of the potential clusters have not been assigned yet
          #   then just do a random one
          fclust <- sample(c(1:ncol(zdat))[aclust==T],1)
        }
      }
      res[val] <- fclust
      #print(res)
    }
  }
  newres <- list(data.frame(cluster=res),
                 data.frame(cluster=c(1:ncol(zdat)),
                            quality=clust.qual(zdat,res)))
  return(newres)
}

##############################################################################
# Cluster Quality is the variance of the the membership similarity for the
#   cluster of those currently (e.g., assigned) belonging to cluster

clust.qual <- function(memb.mat,classign,wtval=NULL){ 
  # memb.mat - results from membership similarity (i.e., memb.sim)
  # classign - cluster assignment for observation can have missing values
  #           for observations that have yet to be assigned to a cluster
  # wtval    - the sampling/IPW weights 
  
  if(is.null(wtval)) wtval <- rep(1,nrow(memb.mat))
  qual <- rep(NA,ncol(memb.mat))
  for(cval in unique(classign[order(classign)])){
    if(!is.na(cval)){
      goodrow <- classign==cval & !is.na(classign)
      qual[cval] <- wtd.var(memb.mat[goodrow,cval],weights=wtval[goodrow])
      #possible that while enforcing hard clustering an observation may be first in cluster
      # in which case a wtd.var with weight 1 returns NA whereas when weight!=1 returns 0.
      # so to allow for hard clustering during assign.clust() will have =0 when NA
      if(is.na(qual[cval])) qual[cval] <- 0
    }
  }
  return(qual)
}

##############################################################################
# Cluster certainty is the mean of the the membership similarity for the 
#   observations that have been assigned to that cluster
clust.certainty <- function(zdat){
  # zdat is membership similarity
  res <- rep(NA,ncol(zdat))
  for(val in c(1:ncol(zdat))){
    res[val] <- mean(zdat[zdat[,val]>0,val])
  }
  return(res)#apply()
}


###########################################################################
# Function to examine cluster merges as alpha is lowered

# Function Input:
  # data -- data set with 1 column per clustering algorithm with cluster assignment for each observation
  # start -- alpha1 value to start with
  # end -- alpha1 value to end with
# Function Output:
  # list with 2 elements
  # n.clusters -- data set with the number of clusters after merging by alpha1
  # merges -- data set with merges that occurred at alpha1

merge.alpha1 <- function(data, start, end) {
  
  new.dat <- data
  
  for (a in (seq(start, end, -0.1))) {
    # Compute cluster similarity matrix
    if (a==start) {
      step1 <- cl.sim(tmp=new.dat, bin.ind=T, wt=NULL)
    } else {
      step1 <- cl.sim(tmp=new.dat, bin.ind=F, wt=rep(1, nrow(data)))
    }
    
    # Determine which clusters merged at that alpha
    similarity <- step1$simularity > a
    if (a==start) {
      n.clusters <- tibble(alpha1=1.0, n_clusters=nrow(step1$simularity))
      merges <- as.data.frame(which(similarity, arr.ind=T))
      if (nrow(merges)>0) {merges$alpha1 <- a}
    } else {
      temp.merges <- as.data.frame(which(similarity, arr.ind=T))
      if (nrow(temp.merges)>0) {temp.merges$alpha1 <- a}
      merges <- bind_rows(merges, temp.merges)
    }
    
    # Merge clusters and save result for next iteration
    step2 <- merge.clust(cls.res=step1, alpha=a)
    new.dat <- step2
    
    # Record number of clusters
    n.clusters.temp <- tibble(alpha1=a, n_clusters=length(step2))
    n.clusters <- bind_rows(n.clusters, n.clusters.temp)
  }
  
  merges <- rename(merges, cluster1=row, cluster2=col)
  ensemble_merges <- list(n.clusters=n.clusters, merges=merges)
  
  return(ensemble_merges)
}


