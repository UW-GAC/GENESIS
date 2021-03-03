

# compute variant-specific inflation-factor for multiple variants
# variants are provided as a matrix, each row provides the allele frequencies 
computeVSIF <- function(eafs, ns, sigmas_sqs){
  
  
  # if eafs is a vector, turn it into a matrix
  if (is.null(dim(eafs))){
    names_eafs <- names(eafs)
    eafs <- matrix(eafs, nrow = 1)
    colnames(eafs) <- names_eafs
  }
  
  k <- ncol(eafs)
  
  # checks
  stopifnot((length(ns) == k ) & (length(sigmas_sqs) ==k) ) 
  stopifnot(all(colnames(eafs) == names(ns)) )
  stopifnot(all(colnames(eafs) == names(sigmas_sqs)))
  
  ## compute sample proportions -- relevant for all variants
  sample_proportion <- ns/sum(ns)
  
  xmat <- data.frame(
    intercept <- rep(1, 3*k),
    geno      <- rep(0:2, k),
    Z         <-  model.matrix( ~factor(rep(1:k, each=3)))[,-1]	
  )
  xmat <- as.matrix(xmat)
  
  beta <- rep(0, k+1)
  pr.z <- sample_proportion[rep(1:k, each=3)]
  
  # repeat for each variant:
  res <- data.frame(SE_true = NA, 
                    SE_naive = NA, 
                    Inflation_factor = rep(NA, nrow(eafs)))
  rownames(res) <- rownames(eafs)
  
  for (i in 1:nrow(eafs)){
    pr.x <- dbinom(xmat[,2], 2, rep(eafs[i,], each=3))
    props <- pr.x*pr.z
    
    Bmat <- t(xmat) %*% diag(props) %*% xmat
    Amat <- t(xmat) %*% diag(props) %*% diag( sigmas_sqs[rep(1:k, each=3)] ) %*% xmat
    
    # sandwich formula for computing the SE of effect size allowing for heterogeneous variances:
    SE_true <- sqrt( (solve(Bmat) %*% Amat %*% solve(Bmat))[2,2] )
    
    # standard SE formula:
    SE_naive <- sqrt( solve(Bmat)[2,2] * sum(props*(sigmas_sqs[rep(1:k, each=3)]) ) )
    
    # return a vector of large-sample true, naive SE, and inflation factor
    res[i,c("SE_true", "SE_naive", "Inflation_factor")] <- 
            c(SE_true, SE_naive, (SE_true/SE_naive)^2)
  }
  
 return(res)
  
}


# a function that gets a null model and extract information from the null model
# group_var_vec is a named vector. Names are sample.ids. Values define groups.
# groups need to correspond to eafs. 
computeVSIFnullmod <- function(nullmod, eafs, group_var_vec){
  
  # if eafs is a vector, turn it into a matrix
  if (is.null(dim(eafs))){
    names_eafs <- names(eafs)
    eafs <- matrix(eafs, nrow = 1)
    colnames(eafs) <- names_eafs
  }
  
  # define groups
  groups <- colnames(eafs)
  
  # check that groups specificied by allele frequencies match
  # the groups specified by group_var_vec
  stopifnot(all(is.element(group_var_vec, groups)))
  
  
  # subset group_var_vec (if needed) to match sample.ids in nullmod 
  if (!is.null(nullmod$sample.id)){
    nullmod_sample.id <- nullmod$sample.id
  } else{ # there is not sample.id entry
    nullmod_sample.id <- rownames(nullmod$model.matrix)
  }
  
  # check that all nullmod_sample.id are in names of group_var_vec
  stopifnot(all(is.element(nullmod_sample.id, names(group_var_vec))))

  # make sure that are ordered appropriately
  group_var_vec <- group_var_vec[nullmod_sample.id]
  
  ## prepare input for function computeVSIF
  ns <- sigmas_sqs <- rep(NA, length(groups))
  names(ns) <- names(sigmas_sqs) <- groups
  
  # extract marginal residuals from the nullmod object:
  if (!is.null(nullmod$resid.marginal)){
    resid <- nullmod$resid.marginal
    } else{
    resid <- nullmod$fit$resid.marginal
  }
  
  for (i in 1:length(groups)){
    ns[[groups[i]]] <- sum(group_var_vec == groups[i])
    sigmas_sqs[[groups[i]]] <- mean(resid[which(group_var_vec == groups[i])]^2)
  }
  
  return(computeVSIF(eafs, ns, sigmas_sqs))
  
  
  
}
