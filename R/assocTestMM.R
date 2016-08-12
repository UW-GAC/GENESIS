assocTestMM <- function(genoData, 
                        nullMMobj,
                        test = "Wald", 
                        snp.include = NULL,
                        chromosome = NULL,
                        impute.geno = TRUE, 
                        snp.block.size = 5000,
                        ivars = NULL,
                        ivar.return.betaCov = FALSE,
                        verbose = TRUE){

    if(class(genoData) == "SeqVarData"){
        # save the filter
        seqFilt.original <- seqGetFilter(genoData)
        # reset so indexing works
        snp.filter <- seqGetData(genoData, "variant.id")
        snp.include <- if(is.null(snp.include)) snp.filter else intersect(snp.include, snp.filter)
        seqResetFilter(genoData, verbose=FALSE)
    }

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
    if(class(genoData) == "GenotypeData"){
        chr <- getChromosome(genoData, index=snp.include$index)
    }else if(class(genoData) == "SeqVarData"){
        chr <- seqGetData(genoData, "chromosome")[snp.include$index]
    }
  
    # X chromosome check for sex variable
    if(.XchromCode(genoData) %in% chr & !.hasSex(genoData)){
        stop("Sex values for the samples are required to compute MAF for X chromosome SNPs")
    }
  
    # Y chromosome 
    if(.YchromCode(genoData) %in% chr){
        # check for sex variable
        if(!.hasSex(genoData)){
            stop("Sex values for the samples are required for Y chromosome SNPs")
        }
        if(class(genoData) == "GenotypeData"){
            getSex(genoData, index = scan.include$index)
        }else if(class(genoData) == "SeqVarData"){
            sex <- sampleData(genoData)$sex[scan.include$index]
        }
        if(!all(sex == "M" )){
            stop("Y chromosome SNPs should be analyzed with only males")
        }
    }

    if(verbose) message("Running analysis with ", scan.include$n, " Samples and ", snp.include$n, " SNPs")

    # GxE interaction
    if(!is.null(ivars)){
        if(test != "Wald"){
            stop("Wald test must be used for Genotype Interaction")
        }
        # names of the interaction variables
        ivnames <- unique(unlist(strsplit(ivars,"[*:]")))
        # index of variables in model.matrix 
        iv.idx <- NULL
        for(val in c("(Intercept)", ivnames)){
            tmp <- grep(val, colnames(nullMMobj$model.matrix))
            if(length(tmp) == 0){
                stop("All of the variables in ivars must have been included as covars in the null model")
            }
            iv.idx <- append(iv.idx, tmp)
        }
        # design matrix for interaction terms
        iW.full <- nullMMobj$model.matrix[,iv.idx]

        # covariance matrix of interaction betas
        if(ivar.return.betaCov){
            res.Vbetas <- vector("list", snp.include$n)
        }
    }  
     
    # set up results matrix
    if(test == "Wald"){
        nv <- c("snpID","chr","n","MAF","minor.allele")
        if(is.null(ivars)){
            nv <- append(nv, c("Est","SE","Wald.Stat","Wald.pval"))
        }else{
            nv <- append(nv, c("Est.G", paste("Est.G",colnames(iW.full)[-1],sep=":"), "SE.G", paste("SE.G",colnames(iW.full)[-1],sep=":"), "GxE.Stat","GxE.pval","Joint.Stat","Joint.pval"))
        }
    }else if(test == "Score"){
        nv <- c("snpID","chr","n","MAF","minor.allele","Score","Var","Score.Stat","Score.pval")
    }
    res <- matrix(NA, nrow=snp.include$n, ncol=length(nv), dimnames=list(NULL, nv))

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
        if(class(genoData) == "GenotypeData"){
            chromChar <- getChromosome(genoData, index = snp.include$index[bidx], char=TRUE)
        }else if(class(genoData) == "SeqVarData"){
            chromChar <- seqGetData(genoData, "chromosome")
        }
        # get sex for selected samples
        if(.hasSex(genoData)){
            if(class(genoData) == "GenotypeData"){
                sex <- getSex(genoData, index = getScanID(genoData) %in% rownames(geno))
            }else if(class(genoData) == "SeqVarData"){
                sex <- sampleData(genoData)$sex
            }
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
        if(b == 1 | !all(keep.previous == keep.geno)){
            if (b > 1) warning("Sample Set Changed: Re-Calculating Matrices!")
      
            # covariate matrix
            W <- nullMMobj$model.matrix[keep.geno, , drop = FALSE]
            k <- ncol(W)
            # outcome
            Y <- nullMMobj$workingY[keep.geno]
            # interaction variables matrix
            if(!is.null(ivars)){
                iW <- iW.full[keep.geno, , drop = FALSE]
                v <- ncol(iW)
            }
      
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
        

        # perform regressions and collect results
        if(test == "Wald"){
            # no interaction
            if(is.null(ivars)){
                Xtilde <- crossprod(Mt,geno) # adjust genotypes for correlation structure and fixed effects
                XtX <- colSums(Xtilde^2) # vector of X^T SigmaInv X (for each SNP)                
                XtX[maf==0] <- NA  # filter monomorphic SNPs
                XtY <- as.vector(crossprod(Xtilde,Ytilde))
                beta <- XtY/XtX
                RSS <- as.numeric( (sY2 - XtY*beta)/(n - k - 1) )
                Vbeta <- RSS/XtX
                Stat <- beta^2/Vbeta

                res[bidx,"Est"] <- beta
                res[bidx,"SE"] <- sqrt(Vbeta)    
                res[bidx,"Wald.Stat"] <- Stat
                res[bidx,"Wald.pval"] <- pchisq(Stat, df=1, lower.tail=FALSE)

            # interaction
            }else{
                GxE.Stat <- rep(NA, ncol(geno))
                Joint.Stat <- rep(NA, ncol(geno))
                # interaction betas and SEs
                res.betas <- matrix(NA, nrow=ncol(geno), ncol=v, dimnames=list(NULL, c("Est.G", paste("Est.G",colnames(iW)[-1],sep=":"))))
                res.SEs <- matrix(NA, nrow=ncol(geno), ncol=v, dimnames=list(NULL, c("SE.G", paste("SE.G",colnames(iW)[-1],sep=":"))))

                for(g in 1:ncol(geno)){
                    # filter monomorphic or missing SNPs
                    if(maf[g] == 0 || is.na(freq[g])){ next }
                    Xtilde <- crossprod(Mt,geno[,g]*iW)
                    XtX <- crossprod(Xtilde)
                    XtXinv <- tryCatch( chol2inv(chol(XtX)), error=function(e){TRUE})  # this is inverse A matrix of sandwich
                    # check that the error function above hasn't been called (which returns TRUE instead of the inverse matrix)
                    if(is.logical(XtXinv)){ next }
                    XtY <- crossprod(Xtilde, Ytilde)
                    betas <- crossprod(XtXinv, XtY)
                    res.betas[g,] <- betas
                    RSS <- as.numeric( (sY2 - crossprod(XtY,betas))/(n - k - v))
                    Vbetas <- XtXinv*RSS
                    if(ivar.return.betaCov){ res.Vbetas[[bidx[g]]] <- Vbetas }  # save covariance matrix for betas?
                    res.SEs[g,] <- sqrt(diag(Vbetas))
                    GxE.Stat[g] <- tryCatch( crossprod(betas[-1],crossprod(chol2inv(chol(Vbetas[-1,-1])),betas[-1])), error=function(e){NA})
                    Joint.Stat[g] <- tryCatch( crossprod(betas,crossprod(XtX,betas))/RSS, error=function(e){NA})
                } # end SNP loop

                res[bidx,c("Est.G", paste("Est.G",colnames(iW)[-1],sep=":"))] <- res.betas
                res[bidx,c("SE.G", paste("SE.G",colnames(iW)[-1],sep=":"))] <- res.SEs
                res[bidx,"GxE.Stat"] <- GxE.Stat
                res[bidx,"GxE.pval"] <- pchisq(GxE.Stat, df=(v-1), lower.tail=FALSE)
                res[bidx,"Joint.Stat"] <- Joint.Stat
                res[bidx,"Joint.pval"] <- pchisq(Joint.Stat, df=v, lower.tail=FALSE)
            } 
      
        }else if(test == "Score"){
            Xtilde <- crossprod(Mt,geno) # adjust genotypes for correlation structure and fixed effects
            XtX <- colSums(Xtilde^2) # vector of X^T P X (for each SNP) b/c (M^T M) = P            
            XtX[maf==0] <- NA  # filter monomorphic SNPs
            score <- as.vector(crossprod(Xtilde,Ytilde)) # X^T P Y
            score[maf==0] <- NA
            Stat <- score^2/XtX

            res[bidx,"Score"] <- score
            res[bidx,"Var"] <- XtX
            res[bidx,"Score.Stat"] <- Stat
            res[bidx,"Score.pval"] <- pchisq(Stat, df=1, lower.tail=FALSE)
        }    
            
        # update keep.previous before moving to next block
        keep.previous <- keep.geno

        # report timing    
        endTime <- Sys.time()
        rate <- format(endTime - startTime, digits=4)    
        if(verbose) message(paste("Block", b, "of", snp.blocks$n, "Completed -", rate))
    } # end block loop
  
    # results data frame
    res <- as.data.frame(res)
  
    # add in snpID
    res[,"snpID"] <- snp.include$value
  
    # chromosome
    res[,"chr"] <- chr

    # convert minor.allele coding back to A/B
    if(class(genoData) == "GenotypeData"){
        res[,"minor.allele"][res[,"minor.allele"] == 1] <- "A"
        res[,"minor.allele"][res[,"minor.allele"] == 0] <- "B"
    }else if(class(genoData) == "SeqVarData"){
        res[,"minor.allele"][res[,"minor.allele"] == 0] <- "ref"
        res[,"minor.allele"][res[,"minor.allele"] == 1] <- "alt"
    }

    # if saving covariance matrix of betas
    if(!is.null(ivars) & ivar.return.betaCov){
        names(res.Vbetas) <- res$snpID
        res <- list(results = res, betaCov = res.Vbetas)
    }
  
    if(class(genoData) == "SeqVarData"){ seqSetFilter(genoData, sample.sel = seqFilt.original$sample.sel, variant.sel = seqFilt.original$variant.sel, verbose = FALSE) }

    return(res)
}

