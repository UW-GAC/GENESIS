pcair <-
function(genoData, v = 10, MAF = 0.05, kinMat = NULL, kin.thresh = 0.025, divMat = NULL, div.thresh = -0.025, unrel.set = NULL, scan.include = NULL, snp.include = NULL, Xchr = FALSE, block.size = 10000, verbose = TRUE){
	
	# checks
	if(MAF < 0 | MAF > 0.5){
		stop("MAF must be between 0 and 0.5")
	}
	
	# load snpIDs
    snpIDs <- getSnpID(genoData)
	# snps to include
	if(!is.null(snp.include)){
		if(!all(snp.include %in% snpIDs)){
			stop("Not all of the SNP IDs in snp.include are in genoData")
		}
    }else{
        if(!Xchr){
            # all autosomal SNPs
            snp.include <- snpIDs[getChromosome(genoData) <= 22]
        }else{
            # all X chromosome SNPs
            snp.include <- snpIDs[getChromosome(genoData) == XchromCode(genoData)]
        }
    }
	# get index of SNP number to include
	snp.include.idx <- which(snpIDs %in% snp.include); rm(snpIDs)
    
    
	# load scanIDs
	scanID <- getScanID(genoData)    
    # samples to include
    if(!is.null(scan.include)){
        if(!all(scan.include %in% scanID)){
            stop("Not all of the scan IDs in scan.include are in genoData")
        }
    }else{
        scan.include <- scanID
    }
    # get index of samples numbers to include
    scan.include.idx <- which(scanID %in% scan.include)
    
    
    # sample size
    nsamp <- length(scan.include)
    if(nsamp == 0){  stop("None of the samples in scan.include are in genoData")  }
      
    # check that scan.include matches kinMat
    if(!is.null(kinMat)){
        if(is.null(colnames(kinMat))){
            stop("colnames and rownames of kinMat must be individual IDs")
        }
        # subset
        kinMat <- subset(kinMat, is.element(rownames(kinMat), scan.include), is.element(colnames(kinMat), scan.include))
        if(!all(scan.include == colnames(kinMat)) | !all(scan.include == rownames(kinMat))){
            stop("colnames and rownames of kinMat must match the scanIDs of genoData")
        }
    }
    
    # check that scan.include matches divMat
    if(!is.null(divMat)){
        if(is.null(colnames(divMat))){
            stop("colnames and rownames of divMat must be individual IDs")
        }
        # subset
        divMat <- subset(divMat, is.element(rownames(divMat), scan.include), is.element(colnames(divMat), scan.include))
        if(!all(scan.include == colnames(divMat)) | !all(scan.include == rownames(divMat))){
            stop("colnames and rownames of divMat must match the scanIDs of genoData")
        }
    }
    
    # check that scan.include matches unrel.set
    if(!is.null(unrel.set)){
        if(!all(unrel.set %in% scan.include)){
            stop("All of the samples in unrel.set must be in the scanIDs of genoData (and scan.include if specified)")
        }
    }
	
	# partition into related and unrelated sets
	if(is.null(kinMat)){
		if(is.null(unrel.set)){
            if(verbose){  message("kinMat and unrel.set both unspecified, Running Standard Principal Components Analysis")  }
			rels <- NULL
			unrels <- scan.include
		}else{
            if(verbose){  message("kinMat not specified, using unrel.set as the Unrelated Set")  }
			rels <- scan.include[!(scan.include %in% unrel.set)]
			unrels <- unrel.set
		}
	}else{
		if(is.null(unrel.set)){
            if(verbose){  message("Partitioning Samples into Related and Unrelated Sets...")  }
		}else{
            if(verbose){  message("Partitioning Samples into Related and Unrelated Sets, unrel.set forced into the Unrelated Set")  }
		}
		part <- pcairPartition(kinMat = kinMat, kin.thresh = kin.thresh, divMat = divMat, div.thresh = div.thresh, unrel.set = unrel.set)
		rels <- part$rels
		unrels <- part$unrels
		if(is.null(rels)){
            if(verbose){  message("No relatives identified, Running Standard Principal Components Analysis")  }
		}
	}
	
	# create index for related and unrealted sets
    # relative to genoData
	unrel.idx <- which(scanID %in% unrels)
    # relative to scan.include
    rel.subidx <- which(scan.include %in% rels)
    unrel.subidx <- which(scan.include %in% unrels)
    # number of samples in each set
	nr <- length(rel.subidx)
	nu <- length(unrel.subidx)
	if(nr > 0){
		method <- "PC-AiR"
	}else{
		method <- "Standard PCA"
	}
    if(verbose){  message(paste("Unrelated Set:",nu,"Samples \nRelated Set:",nr,"Samples"))  }
	
	# determine blocks
	# number of autosomal snps
	nloci <- length(snp.include)
    if(verbose){  message(paste("Running Analysis with",nloci,"SNPs ..."))  }
	# number of blocks of snps
	nblocks <- ceiling(nloci/block.size)
	# start and end positions for SNP blocks
    if(nblocks == 1){
        snp.start <- 1
        snp.end <- nloci
    }else{
        snp.start <- (0:(nblocks-1))*block.size+1
        snp.end <- c( (1:(nblocks-1))*block.size, nloci )
    }
	
	# if relatives
	if(nr > 0){
		# correlation matrix for sPCA
		Psiu <- matrix(0, nrow=nu, ncol=nu)
		# number of snps used
		nsnps <- 0
		
		for(bn in 1:nblocks){
            if(verbose){  message(paste("Computing Genetic Correlation Matrix for the Unrelated Set: Block",bn,"of",nblocks,"..."))  }
			# load genotype data for unrelated set
            geno <- getGenotypeSelection(genoData, snp = snp.include.idx[snp.start[bn]:snp.end[bn]], scan = unrel.idx, drop = FALSE)
            
			# allele freq est from unrelated set
			pA <- 0.5*rowMeans(geno, na.rm = TRUE)
			# remove monomorphic SNPs
			snp.excl <- which(is.na(pA) | pA <= MAF | pA >= (1-MAF))
			if(length(snp.excl > 0)){
				geno <- geno[-snp.excl,]
				pA <- pA[-snp.excl]
			}
			nsnps <- nsnps + length(pA)
			
			# Standardized genotype values
			# estimated variance at each SNP
			sigma.hat <- sqrt(2*pA*(1-pA))
			# z_i = (x_i - 2\hat{p})/sqrt{2*\hat{p}(1-\hat{p})}
			Zu <- (geno-2*pA)/sigma.hat
			Zu[which(is.na(Zu))] <- 0
			
			# unrelated empirical correlation matrix
			Psiu <- Psiu + crossprod(Zu)
		}
		Psiu <- (1/nsnps)*Psiu
		
		# sPCA analysis
        if(verbose){  message("Performing PCA on the Unrelated Set...")  }
		eigu <- eigen(Psiu, symmetric=TRUE)
		
		# subset desired number of eigenvectors
		if(is.null(v)){
			v <- nu
		}
		L <- eigu$values[1:v]
		V <- eigu$vectors[,1:v]
		# sum of eigenvalues
		sum.values <- sum(eigu$values)
		
		# matrix of pseudo-eigenvectors
		Q <- matrix(0, nrow=nr, ncol=v)
		
		# project for related set
        if(verbose){  message("Predicting PC Values for the Related Set...")  }
		for(bn in 1:nblocks){
			# load genotype data
			geno <- getGenotypeSelection(genoData, snp = snp.include.idx[snp.start[bn]:snp.end[bn]], scan = scan.include.idx, drop = FALSE)
            
			# allele freq est from unrelated set
			pA <- 0.5*rowMeans(geno[,unrel.subidx], na.rm = TRUE)
			# remove monomorphic SNPs
			snp.excl <- which(is.na(pA) | pA <= MAF | pA >= (1-MAF))
			if(length(snp.excl > 0)){
				geno <- geno[-snp.excl,]
				pA <- pA[-snp.excl]
			}
			
			# Standardized genotype values
			# estimated variance at each SNP
			sigma.hat <- sqrt(2*pA*(1-pA))
			# z_i = (x_i - 2\hat{p})/sqrt{2*\hat{p}(1-\hat{p})}
			Z <- (geno-2*pA)/sigma.hat
			Z[which(is.na(Z))] <- 0
			
			# subset unrelated and related sets
			Zu <- Z[,unrel.subidx]
			Zr <- Z[,rel.subidx]
			
			# SNP weights
			WT <- tcrossprod(t(V),Zu)
			WL <- (1/L)*WT
			
			# pseudo eigenvectors for relateds
			Q <- Q + crossprod(Zr,t(WL))
		}
		Q <- (1/nsnps)*Q
		
		# concatenate
        if(verbose){  message("Concatenating Results...")  }
		EIG <- matrix(NA, nrow=nsamp, ncol=v)
		EIG[unrel.subidx,] <- V
		EIG[rel.subidx,] <- Q
		
	# if no relatives
	}else{
		# correlation matrix for sPCA
		Psi <- matrix(0, nrow=nsamp, ncol=nsamp)
		# number of snps used
		nsnps <- 0
		
		for(bn in 1:nblocks){
            if(verbose){  message(paste("Computing Genetic Correlation Matrix: Block",bn,"of",nblocks,"..."))  }
			# load genotype data
			geno <- getGenotypeSelection(genoData, snp = snp.include.idx[snp.start[bn]:snp.end[bn]], scan = scan.include.idx, drop = FALSE)

			# allele freq est from entire sample
			pA <- 0.5*rowMeans(geno, na.rm = TRUE)
			# remove monomorphic SNPs
			snp.excl <- which(is.na(pA) | pA <= MAF | pA >= (1-MAF))
			if(length(snp.excl > 0)){
				geno <- geno[-snp.excl,]
				pA <- pA[-snp.excl]
			}
			nsnps <- nsnps + length(pA)
			
			# Standardized genotype values
			# estimated variance at each SNP
			sigma.hat <- sqrt(2*pA*(1-pA))
			# z_i = (x_i - 2\hat{p})/sqrt{2*\hat{p}(1-\hat{p})}
			Z <- (geno-2*pA)/sigma.hat
			Z[which(is.na(Z))] <- 0
			
			# empirical correlation matrix
			Psi <- Psi + crossprod(Z)
		}
		Psi <- (1/nsnps)*Psi
		
		# sPCA analysis
        if(verbose){  message("Performing Standard PCA...")  }
		eig <- eigen(Psi, symmetric=TRUE)
		
		# subset desired number of eigenvectors
		if(is.null(v)){
			v <- nsamp
		}
		# output
		EIG <- eig$vectors[,1:v]
		L <- eig$values[1:v]
		sum.values <- sum(eig$values)
	}
    
    # add scanIDs as rownames of EIG
    rownames(EIG) <- scan.include
	
	# return results
	out <- list(vectors = EIG, 
				values = L, 
				sum.values = sum.values, 
				rels = rels, 
				unrels = unrels,
				kin.thresh = kin.thresh,
				div.thresh = -abs(div.thresh),
				nsamp = nsamp,
				nsnps = nsnps,
				MAF = MAF,
				call = match.call(),
				method = method)
	class(out) <- "pcair"
	return(out)
}
