pcrelate <- function(   genoData, 
                        pcMat = NULL,
                        ibd.probs = TRUE,
                        scan.include = NULL,
                        training.set = NULL,
                        scan.block.size = 5000,
                        snp.include = NULL,
                        Xchr = FALSE,
                        snp.block.size = 10000,
                        MAF = 0.01,
                        write.to.gds = FALSE,
                        gds.prefix = NULL,
                        correct = TRUE,
                        verbose = TRUE){

	# MAF check
	if(MAF < 0 | MAF > 0.5){ stop("MAF must be in [0,0.5]") }

	# SNPs to include in analysis
	snp.include <- getSnpIndex(genoData, snp.include, Xchr)
	# SNP blocks
	snp.blocks <- getBlocks(snp.include$n, snp.block.size)
	if(verbose){ message(paste("Running Analysis with",snp.include$n,"SNPs - in",snp.blocks$n,"Block(s)")) }

	# Scans to include in analysis
	scan.include <- getScanIndex(genoData, scan.include)
    if(scan.include$n == 0){  stop("None of the samples in scan.include are in genoData")  }
	# Scan blocks
    scan.blocks <- getBlocks(scan.include$n, scan.block.size)
    if(verbose){ message(paste("Running Analysis with",scan.include$n,"Samples - in",scan.blocks$n,"Block(s)")) }

    # check number of scan blocks
    if(scan.blocks$n > 1 & !write.to.gds){
    	stop("When the number of Samples is > scan.block.size, output must be written to GDS (i.e. write.to.gds = TRUE) \n For large Samples it is strongly encouraged to write to GDS")
    }

    # Scans to include in training.set
    if(is.null(training.set)){
        # use all scans included in the analysis
        training.set <- scan.include
    }else{
        # subset to training.set
        if(!all(training.set %in% scan.include$value)){
            stop("Not all of the scanIDs in training.set are in scan.include")
        }
        training.set <- getScanIndex(genoData, training.set)
    }


    # PC checks
    if(is.null(pcMat)){
    	if(verbose){ message("pcMat not specified - Calculating Unadjusted Estimates")}
    	method <- "Unadjusted"
        if(training.set$n == scan.include$n){
            if(verbose){ message(paste("Using all",training.set$n,"Samples to Estimate Allele Frequencies")) }
        }else{
            if(verbose){ message(paste("Using",training.set$n,"Samples in training.set to Estimate Allele Frequencies")) }
        }

        # compute part of Hat matrix:
        pcMat <- matrix(1, nrow = scan.include$n, ncol = 1)
        pcProd <- matrix(1/training.set$n, nrow = 1, ncol = training.set$n)

    }else{
        # convert to matrix
        pcMat <- as.matrix(pcMat)
        # check for rownames
        if(is.null(rownames(pcMat))){
            stop("pcMat must have rownames corresponding to ScanIDs")
        }
    	# subset to scan.include
    	pcMat <- pcMat[rownames(pcMat) %in% scan.include$value, , drop=FALSE]
    	if(!all(rownames(pcMat) == scan.include$value)){
    		stop("The order of samples in pcMat does not match the order in genoData")
    	}
    	if(verbose){ message(paste("Using",ncol(pcMat),"PC(s) in pcMat to Calculate Adjusted Estimates")) }
    	method <- "PC-Relate"

        # add intercept
        pcMat <- cbind(rep(1,scan.include$n),pcMat)
        # compute part of the Hat matrix:  (X'X)^{-1}X'
        if(training.set$n == scan.include$n){
            if(verbose){ message(paste("Using all",training.set$n,"Samples to Estimate PC effects on Allele Frequencies")) }
            pcProd <- tcrossprod( chol2inv(chol(crossprod(pcMat))), pcMat )
        }else{
            if(verbose){ message(paste("Using",training.set$n,"Samples in training.set to Estimate PC effects on Allele Frequencies")) }
            pcMat.sub <- pcMat[rownames(pcMat) %in% training.set$value, ]
            pcProd <- tcrossprod( chol2inv(chol(crossprod(pcMat.sub))), pcMat.sub ); rm(pcMat.sub)
        }
    }


    # set up GDS files
    if(write.to.gds){
        # check that gdsfmt is loaded
        requireNamespace("gdsfmt")
    	if(verbose){ message("Creating GDS file for Individual-Specific Allele Frequencies:")}
    	# if no gds.prefix specified
    	if(is.null(gds.prefix)){ gds.prefix <- "tmp" }
    	# create an empty GDS file
    	gds <- createfn.gds(paste0(gds.prefix,"_isaf.gds"))
    	# add sample.id (closed)
    	add.gdsn(gds, "sample.id", scan.include$value, compress="ZIP_RA", closezip=TRUE)
    	# add SNP info (open to append data)
    	add.gdsn(gds, "snp.id", storage="integer", valdim=0, compress="ZIP_RA")
    	add.gdsn(gds, "snp.position", storage="integer", valdim=0, compress="ZIP_RA")
    	add.gdsn(gds, "snp.chromosome", storage="uint8", valdim=0, compress="ZIP_RA")
    	# create node for individual-specific allele frequency data
    	isaf.node <- add.gdsn(gds, "genotype", storage="float32", valdim=c(scan.include$n, 0), compress="ZIP_RA:8M")
    	put.attr.gdsn(isaf.node, "sample.order")
    	sync.gds(gds)

        if(verbose){ message("Estimating Individual-Specific Allele Frequencies...") }
    	# loop through SNP blocks
    	for(bn in 1:snp.blocks$n){    		
    		# keep track of time for rate reporting
    		startTime <- Sys.time()

    		# index of SNPs for this block
    		snp.block.idx <- snp.blocks$start[bn]:snp.blocks$end[bn]

    		# load genotype data
    		geno <- getGenotypeSelection(genoData, snp = snp.include$index[snp.block.idx], scan = training.set$index, drop = FALSE)

    		# estimate individual-specific allele frequencies
    		phat <- estISAF(geno, pcMat, pcProd, transpose = TRUE)

    		# append SNP info
    		append.gdsn(index.gdsn(gds, "snp.id"), val = snp.include$value[snp.block.idx])
    		append.gdsn(index.gdsn(gds, "snp.position"), val = getPosition(genoData, index = snp.include$index[snp.block.idx]))
    		append.gdsn(index.gdsn(gds, "snp.chromosome"), val = getChromosome(genoData, index = snp.include$index[snp.block.idx]))
    		# append individual-specific allele frequency data
    		append.gdsn(index.gdsn(gds, "genotype"), val=phat)
    		sync.gds(gds)

    		endTime <- Sys.time()
    		rate <- format(endTime - startTime, digits=4)
    		if(verbose){ message(paste("...SNP Block",bn,"of",snp.blocks$n,"Completed -", rate)) }
    	}

    	# put gds in readmode
    	readmode.gdsn(index.gdsn(gds, "snp.id"))
    	readmode.gdsn(index.gdsn(gds, "snp.position"))
    	readmode.gdsn(index.gdsn(gds, "snp.chromosome"))
    	readmode.gdsn(isaf.node)        
    	# cleanup
    	filename <- gds$filename
    	closefn.gds(gds)
    	cleanup.gds(filename, verbose=verbose)

    	# read in created GDS
    	isafData <- GenotypeData(GdsGenotypeReader(filename))
    
    	# SNPs to include in the analysis (from isafData)
    	snp.include.isaf <- getSnpIndex(isafData, snp.include$value, Xchr)
    	if(!all.equal(snp.include$value, snp.include.isaf$value)){
    		stop("genoData and isafData are not compatible; perhaps the SNP ordering does not match")
    	}

    	# Scans to include in the analysis (from isafData)
    	scan.include.isaf <- getScanIndex(isafData, scan.include$value)
    	if(!all.equal(scan.include$value, scan.include.isaf$value)){
    		stop("genoData and isafData are not compatible; perhaps the Scan ordering does not match")
    	}


    	if(verbose){ message("Creating GDS file for PC-Relate Results:")}    	
    	# create an empty GDS file for results
    	gds <- createfn.gds(paste0(gds.prefix,"_pcrelate.gds"))
    	# add sample.id (closed)
    	add.gdsn(gds, "sample.id", scan.include$value, storage="uint32", compress="ZIP_RA", closezip=TRUE)
    	# create node for kinship values
    	kinship.node <- add.gdsn(gds, "kinship", storage="float32", valdim=c(scan.include$n, scan.include$n))
    	if(ibd.probs){
    		# create node for k0/k2 values
    		ibd.node <- add.gdsn(gds, "ibd.probs", storage="float32", valdim=c(scan.include$n, scan.include$n))
    	}
    	# create node for number of SNPs used
    	nsnp.node <- add.gdsn(gds, "nsnp", storage="uint32", valdim=c(scan.include$n, scan.include$n))
    	sync.gds(gds)
    }


    # loop through scan blocks
    for(j in 1:scan.blocks$n){
    	# index of scans for block j
    	j.block.idx <- scan.blocks$start[j]:scan.blocks$end[j]
    	j.block.n <- length(j.block.idx)

    	for(i in j:1){
            if(verbose){
                if(scan.blocks$n == 1){ message("Computing PC-Relate Estimates...") }
                else{ message("Computing PC-Relate Estimates for Scan Block Pair (",j,",",i,")...") }
            }
            # index of scans for block i
    		i.block.idx <- scan.blocks$start[i]:scan.blocks$end[i]
    		i.block.n <- length(i.block.idx)

    		# set up empty matrices to hold output
    		kinnum <- matrix(0, nrow=i.block.n, ncol=j.block.n)
    		kindenom <- matrix(0, nrow=i.block.n, ncol=j.block.n)
    		if(ibd.probs){
    			k2num <- matrix(0, nrow=i.block.n, ncol=j.block.n)
    			k2denom <- matrix(0, nrow=i.block.n, ncol=j.block.n)
    			k0num <- matrix(0, nrow=i.block.n, ncol=j.block.n)
    			k0denom <- matrix(0, nrow=i.block.n, ncol=j.block.n)
    		}
    		nsnp <- matrix(0, nrow=j.block.n, ncol=i.block.n)


    		# loop through SNP blocks
    		for(bn in 1:snp.blocks$n){
    			# keep track of time for rate reporting
    			startTime <- Sys.time()

    			# index of SNPs for this block
    			snp.block.idx <- snp.blocks$start[bn]:snp.blocks$end[bn]

    			# load genotype data
                geno <- getGenotypeSelection(genoData, snp = snp.include$index[snp.block.idx], scan = scan.include$index[unique(c(i.block.idx,j.block.idx))], drop = FALSE)
                # index for which columns of geno are in each block                 			
    			i.geno.idx <- which(colnames(geno) %in% scan.include$value[i.block.idx])
    			j.geno.idx <- which(colnames(geno) %in% scan.include$value[j.block.idx])


                # get individual-specific allele frequencies
                if(write.to.gds){
                    # load saved individual specific allele frequency data
                    phat <- getGenotypeSelection(isafData, snp = snp.include.isaf$index[snp.block.idx], scan = scan.include.isaf$index[unique(c(i.block.idx,j.block.idx))], drop = FALSE)
                }else{
                    # estimate individual-specific allele frequencies
                    phat <- estISAF(geno[,scan.include$value %in% training.set$value], pcMat, pcProd, transpose = FALSE)
                }
                # opposite allele freq
                qhat <- 1-phat

                # which snps to filter (missing genotype value or fitted value too big/small)
                filt.idx <- which(is.na(geno) | phat <= MAF | phat >= (1-MAF))
                # set filtered values to 0 (so no contribution)
                geno[filt.idx] <- 0
                phat[filt.idx] <- 0
                qhat[filt.idx] <- 0

    			# product pq (filtered values already 0)
                phatqhat <- phat*qhat
    			

    			# determine the number of snps used for each pair
    			snpcount <- matrix(1, nrow=nrow(geno), ncol=ncol(geno)); snpcount[filt.idx] <- 0
    			nsnp <- nsnp + crossprod(snpcount[,j.geno.idx], snpcount[,i.geno.idx]); rm(snpcount)
    			

    			if(ibd.probs){
    				# create dominance coded matrix	
    				genoD <- matrix(0, nrow = nrow(geno), ncol = ncol(geno))
    				idx0 <- which(geno == 0)
                    genoD[idx0] <- phat[idx0]; rm(idx0)
    				idx2 <- which(geno == 2)    				
    				genoD[idx2] <- qhat[idx2]; rm(idx2)
    				genoD <- genoD - phatqhat
    				# set missing & filtered SNPs to 0 (so no contribution)
    				genoD[filt.idx] <- 0

    				# update k2 numerator
    				k2num <- k2num + crossprod(genoD[,i.geno.idx], genoD[,j.geno.idx]); rm(genoD)
    				# update k2 denominator
    				k2denom <- k2denom + crossprod(phatqhat[,i.geno.idx], phatqhat[,j.geno.idx])

    				# indicator matrices of homozygotes
    				IAA <- geno == 2
    				Iaa <- geno == 0
    				# set missing & filtered SNPs to FALSE (so no contribution)
    				IAA[filt.idx] <- FALSE; Iaa[filt.idx] <- FALSE

    				# update k0 numerator
    				k0num <- k0num + crossprod(IAA[,i.geno.idx], Iaa[,j.geno.idx]) + crossprod(Iaa[,i.geno.idx], IAA[,j.geno.idx])
    				#update k0 denominator
    				k0denom <- k0denom + crossprod(phat[,i.geno.idx]^2, qhat[,j.geno.idx]^2) + crossprod(qhat[,i.geno.idx]^2, phat[,j.geno.idx]^2)
    			}
    			rm(qhat)

    			# update kin denominator
    			kindenom <- kindenom + crossprod(sqrt(phatqhat[,i.geno.idx]), sqrt(phatqhat[,j.geno.idx])); rm(phatqhat)
    			# matrix of regression residuals (filtered values already 0)
    			R <- geno-2*phat; rm(geno); rm(phat)
    			# update kin numerator
    			kinnum <- kinnum + crossprod(R[,i.geno.idx], R[,j.geno.idx]); rm(R)

    			endTime <- Sys.time()
    			rate <- format(endTime - startTime, digits=4)
    			if(verbose){ message(paste("...SNP Block",bn,"of",snp.blocks$n,"Completed -", rate)) }
    		}


    		# compute kinship/inbreeding estimates
    		kin <- kinnum/(4*kindenom)

            # small sample correction
            kincorrect <- NULL
            if(scan.blocks$n == 1 & correct & ncol(pcMat) > 1){
                if(verbose){ message("Performing Small Sample Correction...")}                
                # vectors of kinship and inbreeding values
                kin.tmp <- kin[lower.tri(kin)]
                f.tmp <- 2*diag(kin) - 1
                for(k in 2:ncol(pcMat)){
                    # tmp vector of kinship and inbreeding estimates                                       
                    tmp <- c(kin.tmp, f.tmp)
                    # index of sample to use to find correction
                    idxC <- which(tmp < 2^(-11/2))
                    tmp <- tmp[idxC]

                    # covariance due to ancestry
                    Acov <- tcrossprod(pcMat[,k])
                    Adiag <- diag(Acov)
                    Acov <- Acov[lower.tri(Acov)]                    
                    Avec <- c(Acov, Adiag)[idxC]; rm(idxC)
                    
                    # correction factors
                    vals <- lm(tmp ~ I(Avec))$coef; rm(tmp); rm(Avec)
                    kincorrect <- append(kincorrect, vals)

                    # update
                    kin.tmp <- kin.tmp - vals[1] - vals[2]*Acov; rm(Acov) 
                    f.tmp <- f.tmp - vals[1] - vals[2]*Adiag; rm(Adiag)
                }
                # update results                
                kin[lower.tri(kin)] <- kin.tmp
                kin <- t(kin); kin[lower.tri(kin)] <- kin.tmp
                diag(kin) <- 0.5 + 0.5*f.tmp; rm(f.tmp)
            }

            k2correct <- NULL
            if(ibd.probs){
                # compute k2 estimates
                if(i == j){ 
                    f.j <- 2*diag(kin) - 1
                    fprod <- tcrossprod(f.j)
                }else{
                    f.i <- 2*diag(read.gdsn(kinship.node, start=c(i.block.idx[1], i.block.idx[1]), count=c(i.block.n, i.block.n))) - 1
                    fprod <- tcrossprod(f.i,f.j)
                }
                k2 <- k2num/k2denom - fprod

                # small sample correction
                if(scan.blocks$n == 1 & correct & ncol(pcMat) > 1){
                    # vector of k2 estimates
                    k2.tmp <- k2[lower.tri(k2)]
                    for(k in 2:ncol(pcMat)){
                        # index of sample to use to find correction (unrelateds)
                        idxC <- which(kin.tmp < 2^(-11/2))

                        # covariance due to ancestry
                        Acov <- tcrossprod(pcMat[,k])
                        Acov <- Acov[lower.tri(Acov)]

                        # correction factors
                        vals <- lm(k2.tmp[idxC] ~ I(Acov[idxC]) + I(Acov[idxC]^2))$coef; rm(idxC)
                        k2correct <- append(k2correct, vals)

                        # update
                        k2.tmp <- k2.tmp - vals[1] - vals[2]*Acov - vals[3]*Acov^2; rm(Acov)                        
                    }
                    # index of sample to use to find correction
                    idxC <- which(k2.tmp < 2^(-9/2))

                    # correction factors
                    vals <- lm(k2.tmp[idxC] ~ kin.tmp[idxC])$coef
                    k2correct <- append(k2correct, vals)

                    # update results
                    k2.tmp <- k2.tmp - vals[1] - vals[2]*kin.tmp; rm(kin.tmp)
                    k2[lower.tri(k2)] <- k2.tmp
                    k2 <- t(k2); k2[lower.tri(k2)] <- k2.tmp; rm(k2.tmp)
                }                        

                # compute k0 estimates
                k0 <- 1 - 4*kin + k2
                # index for 1st deg rels
                firstDegidx <- which(kin > 2^(-5/2))
                k0[firstDegidx] <- k0num[firstDegidx]/k0denom[firstDegidx]                
            }

            # set up output
            if(i == j){
                nsnp[upper.tri(nsnp)] <- NA
                if(ibd.probs){
                    k0[lower.tri(k0)] <- k2[lower.tri(k2)]
                    diag(k0) <- NA
                }
            }

            # write output
            if(write.to.gds){
                write.gdsn(kinship.node, kin, start=c(i.block.idx[1], j.block.idx[1]), count=c(i.block.n, j.block.n))
                write.gdsn(nsnp.node, nsnp, start=c(j.block.idx[1], i.block.idx[1]), count=c(j.block.n, i.block.n))
                if(ibd.probs){
                    write.gdsn(ibd.node, k0, start=c(i.block.idx[1], j.block.idx[1]), count=c(i.block.n, j.block.n))
                }
                if(i != j){
                    write.gdsn(kinship.node, t(kin), start=c(j.block.idx[1], i.block.idx[1]), count=c(j.block.n, i.block.n))
                    write.gdsn(nsnp.node, matrix(NA, nrow=i.block.n, ncol=j.block.n), start=c(i.block.idx[1], j.block.idx[1]), count=c(i.block.n, j.block.n))                 
                    if(ibd.probs){
                        write.gdsn(ibd.node, t(k2), start=c(j.block.idx[1], i.block.idx[1]), count=c(j.block.n, i.block.n))
                    }
                }
                sync.gds(gds)
            }       

    	} # i loop
    } # j loop

    if(write.to.gds){
        # add additional information
        add.gdsn(gds, "kincorrect", kincorrect, closezip = TRUE)
        add.gdsn(gds, "k2correct", k2correct, closezip = TRUE)
        add.gdsn(gds, "method", method, storage="character", closezip=TRUE)

        # compress nodes and put in readmode
        filename <- gds$filename
        closefn.gds(gds)
        gds <- openfn.gds(filename, readonly = FALSE)
        compression.gdsn(index.gdsn(gds, "kinship"), compress="ZIP_RA:8M")
        readmode.gdsn(index.gdsn(gds, "kinship"))
        compression.gdsn(index.gdsn(gds, "nsnp"), compress="ZIP_RA:8M")
        readmode.gdsn(index.gdsn(gds, "nsnp"))
    	if(ibd.probs){
    		compression.gdsn(index.gdsn(gds, "ibd.probs"), compress="ZIP_RA:8M")
    		readmode.gdsn(index.gdsn(gds, "ibd.probs"))
    	}        
    	# cleanup    	
    	closefn.gds(gds)
    	close(isafData)
    	cleanup.gds(filename, verbose=verbose)

        out <- paste("PC-Relate output saved to", filename, sep=" ")

    }else{
    	if(ibd.probs){
    		out <- list(sample.id = scan.include$value, kinship = kin, ibd.probs = k0, nsnp = nsnp, kincorrect = kincorrect, k2correct = k2correct, call = match.call(), method = method)
    	}else{
    		out <- list(sample.id = scan.include$value, kinship = kin, nsnp = nsnp, kincorrect = kincorrect, call = match.call(), method = method)
    	}
      class(out) <- "pcrelate"
    }

    return(out)
}



