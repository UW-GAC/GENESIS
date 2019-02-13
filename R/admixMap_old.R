admixMap1 <- function(admixDataList,
                     null.model,
                     #impute.local = TRUE,
                     verbose = TRUE){

  # if admixDataList is GenotypeData (one file), convert to a list
  if(class(admixDataList) == "GenotypeData"){
    admixDataList <- list(admixDataList)
  }

  # how many ancestry components to be tested
  v <- length(admixDataList)

  # if admixDataList doesn't have names, assign them
  if(is.null(names(admixDataList))){
    names(admixDataList) <- paste("Anc",1:v,sep="")
  }

  # get scan.include
  sample.id <- null.model$sample.id
  if (!is.null(sample.id)) {
      sample.index <- match(sample.id, getScanID(admixDataList[[1]]))
  } else {
      sample.index <- match(rownames(null.model$model.matrix),
                            sampleNames(getScanAnnotation(admixDataList[[1]])))
      sample.id <- getScanID(admixDataList[[1]])[sample.index]
  }

  # get chromosome information
  chr <- getChromosome(admixDataList[[1]])

  # check that items match in all elements of admixDataList
  if(v > 1){
    scanID <- getScanID(admixDataList[[1]])
    for(i in 2:v){
      # scanIDs
      if(!all(scanID == getScanID(admixDataList[[i]]) )){
        stop("scanIDs do not match for all elements of admixDataList")
      }
      # chromosomes
      if(!all(chr == getChromosome(admixDataList[[i]]) )){
        stop("Chromosomes do not match for all elements of admixDataList")
      }
    }
  }

  # Sex chromosome checks
  Xcode <- rep(NA,v); Ycode <- rep(NA,v)
  for(i in 1:v){
    Xcode[i] <- XchromCode(admixDataList[[i]])
    Ycode[i] <- YchromCode(admixDataList[[i]])
  }
  if(any(Xcode %in% chr)){
    # check for matching
    if(length(unique(Xcode)) != 1){
      stop("X chromosome codes do not match for all elements of admixDataList")
    }
    # check for sex variable
    for(i in 1:v){
      if(!hasSex(admixDataList[[i]])){
        stop("Sex values for the samples are required for X chromosome SNPs")
      }
    }
  }
  if(any(Ycode %in% chr)){
    # check for matching
    if(length(unique(Ycode)) != 1){
      stop("Y chromosome codes do not match for all elements of admixDataList")
    }
    if(!all(chr == Ycode[1])){
      stop("Y chromosome must be analyzed separately")
    }
    # check for sex variable
    for(i in 1:v){
      if(!hasSex(admixDataList[[i]])){
        stop("Sex values for the samples are required for Y chromosome SNPs")
      }
      # check for matching
      if(!all( getSex(admixDataList[[1]]) == getSex(admixDataList[[i]]) )){
        stop("Sex values for the samples do not match for all elements of admixDataList (required for Y chr)")
      }
    }    
    # check for males
    if (!all(getSex(admixDataList[[1]], index=sample.index) == "M")){
        stop("Y chromosome SNPs should be analyzed with only males")
    }
  }  

  
  #if(verbose) message("Running analysis with ", scan.include$n, " Samples and ", snp.include$n, " SNPs")

  
  # set up results matrix
  # add in frequency of each ancestry at the SNP
  res.list <- list()
  nv <- c("snpID","chr","n")
  if(v > 1){
    for(i in 1:v){
      nv <- append(nv, paste(names(admixDataList)[i],c(".freq",".Est",".SE"), sep=""))
    }
    nv <- append(nv, c("Joint.Stat", "Joint.pval"))
  }else{
    nv <- append(nv, c("freq", "Est", "SE", "Stat", "pval"))
  }

  
  n.iter <- length(snpFilter(admixDataList[[1]]))
  b <- 1
  iterate <- TRUE
  while (iterate) {

  res <- matrix(NA, nrow=length(chr), ncol=length(nv), dimnames=list(NULL, nv))
    
  #if(verbose) message("Beginning Calculations...")
    
    # get local ancestry for the block
    local <- array(NA, dim=c(length(sample.id),length(chr),v)) # indices: scan, snp, ancestry
    for(i in 1:v){
      local[,,i] <- getGenotypeSelection(admixDataList[[i]], scanID=sample.id, order="selection", transpose=TRUE)
    }

    # ancestral frequency
    freq <- 0.5*colMeans(local, dims=1, na.rm=T) # matrix:  rows are SNPs, columns are ancestries
    # for X chr
    chr <- getChromosome(admixDataList[[1]])
    if(XchromCode(admixDataList[[1]]) %in% chr){
      # which are on X chr
      Xidx <- chr == XchromCode(admixDataList[[1]])
      # males
      m <- (getSex(admixDataList[[1]]) == "M")
      f <- (getSex(admixDataList[[1]]) == "F")
      # calculate allele freq for X
      freq[Xidx,] <- (0.5*colSums(local[m,Xidx,], dims=1, na.rm=T) + colSums(local[f,Xidx,], dims=1, na.rm=T)) / (colSums(!is.na(local[m,Xidx,]), dims=1) + 2*colSums(!is.na(local[f,Xidx,]), dims=1) )
    }
    
    # add to output
    if(v > 1){
      for(i in 1:v){
        res[,paste(names(admixDataList)[i],".freq", sep="")] <- freq[,i]
      }
    }else{
      freq <- as.vector(freq)
      res[,"freq"] <- freq
    }
    
    # impute missing local ancestry values
    #if(impute.local){      
    #  miss.idx <- which(is.na(local))
    #  if(length(miss.idx) > 0){
    #    freq.idx <- ceiling(miss.idx/n)
    #    local[miss.idx] <- 2*freq[freq.idx]  # double check that this indexes in the array correctly (pretty sure it does)
    #  }
    #}
    # this is imputing as population average.  would it be better to impute to individual's global ancestry?
    ### no missing data currently - worry about this later ###
    
    # sample size
    n <- length(sample.id)
    res[, "n"] <- n
    
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
  res <- as.data.frame(res)
  
  # add in snpID
  res$snpID <- getSnpID(admixDataList[[1]])
    
  # chromosome
  res[,"chr"] <- chr
  
    res.list[[b]] <- res
    
    if (verbose & b %% 100 == 0) {
        message(paste("Iteration", b , "of", n.iter, "completed"))
    }
    b <- b + 1
    for (i in 1:v) {
        iterate <- GWASTools::iterateFilter(admixDataList[[i]])
    }
  }
  
  do.call(rbind, res.list)
}
