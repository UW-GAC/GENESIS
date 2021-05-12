

# compute variant-specific inflation-factor for multiple variants
# variants are provided as a matrix, each row provides the allele frequencies 
computeVSIF <- function(freq, n, sigma.sq){
  
  # if freq is a vector, turn it into a matrix
  if (is.null(dim(freq))){
    names_freq <- names(freq)
    freq <- matrix(freq, nrow = 1)
    colnames(freq) <- names_freq
  }
  
  k <- ncol(freq)
  
  # checks
  stopifnot((length(n) == k ) & (length(sigma.sq) ==k) ) 
  stopifnot(all(colnames(freq) == names(n)) )
  stopifnot(all(colnames(freq) == names(sigma.sq)))
  
  ## compute sample proportions -- relevant for all variants
  sample_proportion <- n/sum(n)
  
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
                    Inflation_factor = rep(NA, nrow(freq)))
  rownames(res) <- rownames(freq)
  
  for (i in 1:nrow(freq)){
    pr.x <- dbinom(xmat[,2], 2, rep(freq[i,], each=3))
    props <- pr.x*pr.z
    
    Bmat <- t(xmat) %*% diag(props) %*% xmat
    Amat <- t(xmat) %*% diag(props) %*% diag( sigma.sq[rep(1:k, each=3)] ) %*% xmat
    
    # sandwich formula for computing the SE of effect size allowing for heterogeneous variances:
    SE_true <- sqrt( (solve(Bmat) %*% Amat %*% solve(Bmat))[2,2] )
    
    # standard SE formula:
    SE_naive <- sqrt( solve(Bmat)[2,2] * sum(props*(sigma.sq[rep(1:k, each=3)]) ) )
    
    # return a vector of large-sample true, naive SE, and inflation factor
    res[i,c("SE_true", "SE_naive", "Inflation_factor")] <- 
            c(SE_true, SE_naive, (SE_true/SE_naive)^2)
  }
  
 return(res)
  
}


# a function that gets a null model and extract information from the null model
# group.var.vec is a named vector. Names are sample.ids. Values define groups.
# groups need to correspond to freq. 
computeVSIFNullModel <- function(null.model, freq, group.var.vec){
    
  # Convert old null model format if necessary.
  null.model <- .updateNullModelFormat(null.model)
  
  # if freq is a vector, turn it into a matrix
  if (is.null(dim(freq))){
    names_freq <- names(freq)
    freq <- matrix(freq, nrow = 1)
    colnames(freq) <- names_freq
  }
  
  # define groups
  groups <- colnames(freq)
  
  # check that groups specificied by allele frequencies match
  # the groups specified by group.var.vec
  stopifnot(all(is.element(group.var.vec, groups)))
  
  
  # subset group.var.vec (if needed) to match sample.ids in null.model 
  if (!is.null(null.model$sample.id)){
    nullmod_sample.id <- null.model$fit$sample.id
  } else{ # there is not sample.id entry
    nullmod_sample.id <- rownames(null.model$model.matrix)
  }
  
  # check that all nullmod_sample.id are in names of group.var.vec
  stopifnot(all(is.element(nullmod_sample.id, names(group.var.vec))))

  # make sure that are ordered appropriately
  group.var.vec <- group.var.vec[nullmod_sample.id]
  
  ## prepare input for function computeVSIF
  n <- sigma.sq <- rep(NA, length(groups))
  names(n) <- names(sigma.sq) <- groups
  
  # extract marginal residuals from the null.model object:
  resid <- null.model$fit$resid.marginal
  
  for (i in 1:length(groups)){
    n[[groups[i]]] <- sum(group.var.vec == groups[i])
    sigma.sq[[groups[i]]] <- mean(resid[which(group.var.vec == groups[i])]^2)
  }
  
  return(computeVSIF(freq, n, sigma.sq))
  
}
