pcair <- function(	genoData, 
					v = 20,					 
					kinMat = NULL, 
					kin.thresh = 2^(-11/2), 
					divMat = NULL, 
					div.thresh = -2^(-11/2), 
					unrel.set = NULL, 
					scan.include = NULL, 
					snp.include = NULL, 
					chromosome = NULL,
					snp.block.size = 10000, 
					MAF = 0.01,
					verbose = TRUE){
	
	# MAF check
	if(MAF < 0 | MAF > 0.5){ stop("MAF must be in [0,0.5]") }

	# SNPs to include in analysis
	snp.include <- getSnpIndex(genoData, snp.include, chromosome)
	# SNP blocks
	snp.blocks <- getBlocks(snp.include$n, snp.block.size)

	# Scans to include in analysis
	scan.include <- getScanIndex(genoData, scan.include)
	if(scan.include$n == 0){  stop("None of the samples in scan.include are in genoData")  }
      
    # check that scan.include matches kinMat
    if(!is.null(kinMat)){
        if(is.null(colnames(kinMat))){
            stop("colnames and rownames of kinMat must be individual IDs")
        }
        # subset
        kinMat <- subset(kinMat, rownames(kinMat) %in% scan.include$value, colnames(kinMat) %in% scan.include$value)
        if(!all(scan.include$value == colnames(kinMat)) | !all(scan.include$value == rownames(kinMat))){
            stop("colnames and rownames of kinMat must match the scanIDs of genoData")
        }
    }
    
    # check that scan.include matches divMat
    if(!is.null(divMat)){
        if(is.null(colnames(divMat))){
            stop("colnames and rownames of divMat must be individual IDs")
        }
        # subset
        divMat <- subset(divMat, rownames(divMat) %in% scan.include$value, colnames(divMat) %in% scan.include$value)
        if(!all(scan.include$value == colnames(divMat)) | !all(scan.include$value == rownames(divMat))){
            stop("colnames and rownames of divMat must match the scanIDs of genoData")
        }
    }
    
    # check that scan.include matches unrel.set
    if(!is.null(unrel.set)){
        if(!all(unrel.set %in% scan.include$value)){
            stop("All of the scanIDs in unrel.set must be in the scanIDs of genoData (and scan.include if specified)")
        }
    }
	
	# partition into related and unrelated sets
	if(is.null(kinMat)){
		if(is.null(unrel.set)){
            if(verbose){  message("kinMat and unrel.set both unspecified, Running Standard Principal Components Analysis")  }
			rels <- NULL
			unrels <- scan.include$value
		}else{
            if(verbose){  message("kinMat not specified, using unrel.set as the Unrelated Set")  }
			rels <- scan.include$value[!(scan.include$value %in% unrel.set)]
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
	unrel.idx <- which(getScanID(genoData) %in% unrels)
    # relative to scan.include
    rel.subidx <- which(scan.include$value %in% rels)
    unrel.subidx <- which(scan.include$value %in% unrels)
    # number of samples in each set
	nr <- length(rel.subidx)
	nu <- length(unrel.subidx)
	if(nr > 0){
		method <- "PC-AiR"
	}else{
		method <- "Standard PCA"
	}
    if(verbose){  message(paste("Unrelated Set:",nu,"Samples \nRelated Set:",nr,"Samples"))  }	
	

	if(verbose){ message(paste("Running Analysis with",snp.include$n,"SNPs - in",snp.blocks$n,"Block(s)")) }
	# if relatives
	if(nr > 0){
		# correlation matrix for sPCA
		Psiu <- matrix(0, nrow=nu, ncol=nu)
		# number of snps used
		nsnps <- 0
		
		for(bn in 1:snp.blocks$n){
            if(verbose){  message(paste("Computing Genetic Correlation Matrix for the Unrelated Set: Block",bn,"of",snp.blocks$n,"..."))  }

            # index of SNPs for this block
            snp.block.idx <- snp.blocks$start[bn]:snp.blocks$end[bn]

			# load genotype data for unrelated set
            geno <- getGenotypeSelection(genoData, snp = snp.include$index[snp.block.idx], scan = unrel.idx, drop = FALSE)
            
			# allele freq est from unrelated set
			pA <- 0.5*rowMeans(geno, na.rm = TRUE)
			# filter SNPs on MAF
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
			# impute to mean
			Zu[which(is.na(Zu))] <- 0
			
			# unrelated empirical correlation matrix
			Psiu <- Psiu + crossprod(Zu)
		}
		Psiu <- (1/nsnps)*Psiu
		
		# PCA
        if(verbose){  message("Performing PCA on the Unrelated Set...")  }
		eigu <- eigen(Psiu, symmetric=TRUE)
		
		# subset desired number of eigenvectors
		if(is.null(v)){ v <- nu	}
		L <- eigu$values[1:v]
		V <- eigu$vectors[,1:v]
		# sum of eigenvalues
		sum.values <- sum(eigu$values)
		
		# matrix of pseudo-eigenvectors
		Q <- matrix(0, nrow=nr, ncol=v)
		
		# project for related set        
		for(bn in 1:snp.blocks$n){
			if(verbose){  message(paste("Predicting PC Values for the Related Set: Block",bn,"of",snp.blocks$n,"..."))  }

			# index of SNPs for this block
            snp.block.idx <- snp.blocks$start[bn]:snp.blocks$end[bn]

			# load genotype data
			geno <- getGenotypeSelection(genoData, snp = snp.include$index[snp.block.idx], scan = scan.include$index, drop = FALSE)
            
			# allele freq est from unrelated set
			pA <- 0.5*rowMeans(geno[,unrel.subidx], na.rm = TRUE)
			# filter SNPs on MAF
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
			# impute to mean
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
		EIG <- matrix(NA, nrow=scan.include$n, ncol=v)
		EIG[unrel.subidx,] <- V
		EIG[rel.subidx,] <- Q
		
	# if no relatives
	}else{
		# correlation matrix for sPCA
		Psi <- matrix(0, nrow=scan.include$n, ncol=scan.include$n)
		# number of snps used
		nsnps <- 0
		
		for(bn in 1:snp.blocks$n){
            if(verbose){  message(paste("Computing Genetic Correlation Matrix: Block",bn,"of",snp.blocks$n,"..."))  }

            # index of SNPs for this block
            snp.block.idx <- snp.blocks$start[bn]:snp.blocks$end[bn]

			# load genotype data
			geno <- getGenotypeSelection(genoData, snp = snp.include$index[snp.block.idx], scan = scan.include$index, drop = FALSE)

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
		if(is.null(v)){	v <- scan.include$n }
		# output
		EIG <- eig$vectors[,1:v]
		L <- eig$values[1:v]
		sum.values <- sum(eig$values)
	}
    
    # add scanIDs as rownames of EIG
    rownames(EIG) <- scan.include$value
	
	# return results
	out <- list(vectors = EIG, 
				values = L, 
				sum.values = sum.values, 
				rels = rels, 
				unrels = unrels,
				kin.thresh = kin.thresh,
				div.thresh = -abs(div.thresh),
				nsamp = scan.include$n,
				nsnps = nsnps,
				MAF = MAF,
				call = match.call(),
				method = method)
	class(out) <- "pcair"
	return(out)
}
