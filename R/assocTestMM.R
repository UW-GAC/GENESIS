assocTestMM <- function(genoData, 
                        nullMMobj, 
                        test = "Wald", 
                        snp.include = NULL,
                        chromosome = NULL,
                        impute.geno = TRUE, 
                        snp.block.size = 5000, 
                        verbose = TRUE){
  
  # check that test is valid
  if(!is.element(test,c("Wald","Score"))){
    stop("test must be one of Wald or Score")
  }

  if(test == "Wald" & nullMMobj$family$family != "gaussian"){
    stop("Wald test can only be performed when the family of the null model is gaussian")
  }

  # get scan.include
  scan.include <- getScanIndex(data = genoData, nullMMobj$scanID)

  # set SNPs to include
  snp.include <- getSnpIndex(data = genoData, snp.include, chromosome)  
  # get chromosome 
  chr <- getChromosome(genoData, index=snp.include$index)
  
  # X chromosome check for sex variable
  if(XchromCode(genoData) %in% chr & !hasSex(genoData)){
    stop("Sex values for the samples are required to compute MAF for X chromosome SNPs")
  }
  
  # Y chromosome 
  if(YchromCode(genoData) %in% chr){
    # check for sex variable
    if(!hasSex(genoData)){
      stop("Sex values for the samples are required for Y chromosome SNPs")
    }
    if(!all( getSex(genoData, index = scan.include$index) == "M" )){
      stop("Y chromosome SNPs should be analyzed with only males")
    }
  }  
  
  if(verbose) message("Running analysis with ", scan.include$n, " Samples and ", snp.include$n, " SNPs")

   
  # set up results matrix
  if(test == "Wald"){
    nv <- c("snpID","chr","n","MAF","minor.allele","Est","SE","Wald.Stat","Wald.pval")
  }else if(test == "Score"){
    nv <- c("snpID","chr","n","MAF","minor.allele","Score","Var","Score.Stat","Score.pval")
  }
  res <- matrix(NA, nrow=snp.include$n, ncol=length(nv), dimnames=list(NULL, nv))
    
  # chromosome
  res[,"chr"] <- chr

  # SNP blocks
  snp.blocks <- getBlocks(snp.include$n, snp.block.size)

  # since we haven't done any loops here yet:
  keep.previous <- rep(TRUE, scan.include$n)

  
  if(verbose) message("Beginning Calculations...")
  # loop through blocks
  for(b in 1:snp.blocks$n){
    
    # keep track of time for rate reporting
    startTime <- Sys.time()

    # set index for this block of SNPs
    bidx <- snp.blocks$start[b]:snp.blocks$end[b]   

    # prepare the genotype data (read in genotypes; mean impute missing values or exclude samples with complete missingness) 
    geno <- prepareGenotype(genoData = genoData, snp.read.idx = snp.include$index[bidx], scan.read.idx = scan.include$index, impute.geno = impute.geno)
    # rows are samples, columns are SNPs

    # samples kept for this block
    keep.geno <- scan.include$value %in% rownames(geno)
    # sample size for this block
    n <- sum(keep.geno)
    res[bidx, "n"] <- n

    # get chromosome for selected SNPs
    chromChar <- getChromosome(genoData, index = snp.include$index[bidx], char=TRUE)
    # get sex for selected samples
    if(hasSex(genoData)){
      sex <- getSex(genoData, index = getScanID(genoData) %in% rownames(geno))
    }else{
      sex <- NULL
    }
    # calculate allele frequencies
    freq <- alleleFreq(geno = geno, chromChar = chromChar, sex = sex)
    # minor allele frequency
    maf <- ifelse(freq < 0.5, freq, 1-freq)
    res[bidx,"MAF"] <- maf
    # minor allele coding:  A = 1, B = 0
    res[bidx,"minor.allele"] <- ifelse(freq < 0.5, 1, 0)


    # matrices for calculating tests
    if (b == 1 | !all(keep.previous == keep.geno)){      
      if (b > 1) warning("Sample Set Changed: Re-Calculating Matrices!")
      
      # covariate matrix
      W <- nullMMobj$model.matrix[keep.geno, , drop=FALSE]
      k <- ncol(W)
      # outcome
      Y <- nullMMobj$workingY[keep.geno]
      
      # here we have to subset the matrix, not the inverse
      # this is a fancy way of getting the inverse of the subset without having to get the original matrix
      # cholesky decomposition of sigma inverse (inverse phenotype covariance matrix)
      chol.idx <- which(!(colnames(nullMMobj$cholSigmaInv) %in% rownames(geno)))
      C <- subsetCholSigmaInv(nullMMobj$cholSigmaInv, chol.idx)

      CW <- crossprod(C, W)
      # matrix used to adjust phenotype and genotype for fixed effect covariates AND "decorrelating" the phenotype and genotype 
      Mt <- C - tcrossprod(tcrossprod(C,tcrossprod(chol2inv(chol(crossprod(CW))),CW)),CW)
      # phenotype adjusted for the covariates/correlation structure
      Ytilde <- crossprod(Mt,Y)
      sY2 <- sum(Ytilde^2)
    }    
        
    # perform regressions      
    if(test == "Wald"){
      Xtilde <- crossprod(Mt,geno) # adjust genotypes for correlation structure and fixed effects
      XtX <- colSums(Xtilde^2) # vector of X^T SigmaInv X (for each SNP)
      # filter monomorphic SNPs
      XtX[maf==0] <- NA
      XtY <- as.vector(crossprod(Xtilde,Ytilde))
      beta <- XtY/XtX
      RSS <- as.numeric( (sY2 - XtY*beta)/(n - k - 1) )
      Vbeta <- RSS/XtX
      Stat <- beta^2/Vbeta
      
    }else if(test == "Score"){
      Xtilde <- crossprod(Mt,geno) # adjust genotypes for correlation structure and fixed effects
      XtX <- colSums(Xtilde^2) # vector of X^T P X (for each SNP) b/c (M^T M) = P
      # filter monomorphic SNPs
      XtX[maf==0] <- NA
      score <- as.vector(crossprod(Xtilde,Ytilde)) # X^T P Y
      score[maf==0] <- NA
      Stat <- score^2/XtX

    }
    
    # collect results
    if(test == "Wald"){
      res[bidx,"Est"] <- beta
      res[bidx,"SE"] <- sqrt(Vbeta)
      #res[bidx,"RSS"] <- RSS
      res[bidx,"Wald.Stat"] <- Stat
      res[bidx,"Wald.pval"] <- pchisq(Stat, df=1, lower.tail=FALSE)
      
    }else if(test == "Score"){
      res[bidx,"Score"] <- score
      res[bidx,"Var"] <- XtX
      res[bidx,"Score.Stat"] <- Stat
      res[bidx,"Score.pval"] <- pchisq(Stat, df=1, lower.tail=FALSE)
    }
    
    endTime <- Sys.time()
    rate <- format(endTime - startTime, digits=4)
    
    # update keep.previous before moving to next block
    keep.previous <- keep.geno
    
    if(verbose) message(paste("Block", b, "of", snp.blocks$n, "Completed -", rate))
  } # end block loop
  
  # results data frame
  res <- as.data.frame(res)
  
  # add in snpID
  res[,"snpID"] <- snp.include$value
  
  # convert minor.allele coding back to A/B
  res[,"minor.allele"][res[,"minor.allele"] == 1] <- "A"
  res[,"minor.allele"][res[,"minor.allele"] == 0] <- "B"
  
  return(res)
}

