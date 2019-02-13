admixMap <- function(admixDataList,
                     null.model,
                     verbose = TRUE){

  # if admixDataList is one file, convert to a list
  if(!is.list(admixDataList)){
    admixDataList <- list(admixDataList)
  }

  # how many ancestry components to be tested
  v <- length(admixDataList)

  # if admixDataList doesn't have names, assign them
  if(is.null(names(admixDataList))){
    names(admixDataList) <- paste("Anc",1:v,sep="")
  }

  # get sample index
    if (is(admixDataList[[1]], "GenotypeIterator")) {
        sample.index <- lapply(admixDataList, .sampleIndexNullModel, null.model)
        if (!.listIdentical(sample.index)) stop("sample IDs do not match for all elements of admixDataList")
        sample.index <- sample.index[[1]]
        sample.id <- names(sample.index)
    } else if (is(admixDataList[[1]], "SeqVarIterator")) {
        sample.index <- lapply(admixDataList, .setFilterNullModel, null.model, verbose=verbose)
        if (!.listIdentical(sample.index)) stop("sample IDs do not match for all elements of admixDataList")
        sample.index <- sample.index[[1]]
    } else {
        stop("admixDataList must contain GenotypeIterator or SeqVarIterator objects")
    }
    n.samp <- length(sample.index)

  # get variant information
    var.info <- lapply(admixDataList, variantInfo)
    if (!.listIdentical(var.info)) stop("variants do not match for all elements of admixDataList")
    var.info <- var.info[[1]]
    n.var <- nrow(var.info)
               
  
  # set up results matrix
  # add in frequency of each ancestry at the SNP
  res.list <- list()
  nv <- c("n.obs")
  if(v > 1){
    for(i in 1:v){
      nv <- append(nv, paste(names(admixDataList)[i],c(".freq",".Est",".SE"), sep=""))
    }
    nv <- append(nv, c("Joint.Stat", "Joint.pval"))
  }else{
    nv <- append(nv, c("freq", "Est", "SE", "Stat", "pval"))
  }

  
  n.iter <- length(variantFilter(admixDataList[[1]]))
  b <- 1
  iterate <- TRUE
  while (iterate) {

  res <- matrix(NA, nrow=n.var, ncol=length(nv), dimnames=list(NULL, nv))
    
    # get local ancestry for the block
    local <- array(NA, dim=c(n.samp,n.var,v)) # indices: scan, snp, ancestry
      for(i in 1:v){
          if (is(admixDataList[[1]], "GenotypeIterator")) {
              local[,,i] <- getGenotypeSelection(admixDataList[[i]], scanID=sample.id, order="selection", transpose=TRUE, use.names=FALSE, drop=FALSE)
          } else {
              local[,,i] <- refDosage(admixDataList[[i]], use.names=FALSE)[sample.index,,drop=FALSE]
          }
      }

    # ancestral frequency
    if(v > 1){
      for(i in 1:v){
          freq <- .alleleFreq(admixDataList[[i]], local[,,i], sample.index=sample.index)
          col <- if (v > 1) paste(names(admixDataList)[i],".freq", sep="") else "freq"
          res[,col] <- freq
      }
    }
    
    # sample size
    res[, "n.obs"] <- n.samp
    
    k <- ncol(null.model$model.matrix)
    Ytilde <- null.model$Ytilde
    sY2 <- sum(Ytilde^2)
        
    # perform regressions
    if(v == 1){
      local <- local[,,1]
      Gtilde <- calcGtilde(null.model, local)
      GPG <- colSums(Gtilde^2) # vector of G^T P G (for each SNP)
      # filter monomorphic SNPs
      GPG[which(freq == 0 || freq == 1)] <- NA
      beta <- as.vector(crossprod(Gtilde,Ytilde)/GPG)
      Vbeta <- (sY2/GPG - beta^2)/(n - k - 1) # RSS/GPG
      Stat <- beta^2/Vbeta

      res[,"Est"] <- beta
      res[,"SE"] <- sqrt(Vbeta)
      res[,"Stat"] <- Stat
      res[,"pval"] <- pchisq(Stat, df=1, lower.tail=FALSE)      

    }else{
      Joint.Stat <- rep(NA, ncol(local))
      Est <- matrix(NA, nrow=ncol(local), ncol=v)
      SE <- matrix(NA, nrow=ncol(local), ncol=v)
      for(g in 1:ncol(local)){
        if(identical(local[,g,], local[,(g-1),])){
          Joint.Stat[g] <- Joint.Stat[g-1]
          Est[g,] <- Est[g-1,]
          SE[g,] <- SE[g-1,]                     
          next
        }
        # filter monomorphic or missing SNPs
        if(any(freq[g,]==1) || sum(freq[g,]==0)){ next }
        Gtilde <- calcGtilde(null.model, local[,g,])
        Ytilde <- null.model$Ytilde
        sY2 <- sum(Ytilde^2)
        GPG <- crossprod(Gtilde)
        GPGinv <- tryCatch( chol2inv(chol(GPG)), error=function(e){TRUE})
        # check that the error function above hasn't been called (which returns TRUE instead of the inverse matrix)
        if(is.logical(GPGinv)){ next }
        GPY <- crossprod(Gtilde, Ytilde)
        betas <- crossprod(GPGinv, GPY)
        RSS <- as.numeric((sY2 - crossprod(GPY,betas))/(n - k - v))
        Vbetas <- GPGinv*RSS

        Est[g,] <- as.vector(betas)
        SE[g,] <- as.vector(sqrt(diag(Vbetas)))
        Joint.Stat[g] <- tryCatch( crossprod(betas,crossprod(GPG,betas))/RSS, error=function(e){NA})
      } # g loop

      # collect results
      for(i in 1:v){
        res[,paste(names(admixDataList)[i],".Est", sep="")] <- Est[,i]
        res[,paste(names(admixDataList)[i],".SE", sep="")] <- SE[,i]
      }
      res[,"Joint.Stat"] <- Joint.Stat
      res[,"Joint.pval"] <- pchisq(Joint.Stat, df=v, lower.tail=FALSE)
    } # else

  # results data frame
  res <- cbind(var.info, res)
  
    res.list[[b]] <- res
    
    if (verbose & b %% 100 == 0) {
        message(paste("Iteration", b , "of", n.iter, "completed"))
    }
    b <- b + 1
    for (i in 1:v) {
        if (is(admixDataList[[i]], "GenotypeIterator")) {
            iterate <- GWASTools::iterateFilter(admixDataList[[i]])
        } else {
            iterate <- SeqVarTools::iterateFilter(admixDataList[[i]])
        }
    }
  }
  
  do.call(rbind, res.list)
}
