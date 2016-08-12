pcrelate <- function(   genoData, 
                        pcMat = NULL,
                        freq.type = "individual",
                        scale = "overall",
                        ibd.probs = TRUE,
                        scan.include = NULL,
                        training.set = NULL,
                        scan.block.size = 5000,
                        snp.include = NULL,
                        chromosome = NULL,
                        snp.block.size = 10000,
                        MAF = 0.01,
                        write.to.gds = FALSE,
                        gds.prefix = NULL,
                        correct = TRUE,
                        verbose = TRUE){

    if(class(genoData) == "SeqVarData"){
        # save the filter
        seqFilt.original <- seqGetFilter(genoData)
        # reset so indexing works
        snp.filter <- seqGetData(genoData, "variant.id")
        snp.include <- if(is.null(snp.include)) snp.filter else intersect(snp.include, snp.filter)
        scan.filter <- seqGetData(genoData, "sample.id")
        scan.include <- if(is.null(scan.include)) scan.filter else intersect(scan.include, scan.filter)
        seqResetFilter(genoData, verbose=FALSE)
    }

    # MAF check
    if(MAF < 0 | MAF > 0.5){ stop("MAF must be in [0,0.5]") }

    # SNPs to include in analysis
    snp.include <- getSnpIndex(genoData, snp.include, chromosome)
    if(snp.include$n == 0){  stop("None of the SNPs in snp.include are in genoData")  }
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

    if(scale == "none" & ibd.probs == TRUE){
        stop("IBD probabilities can only be calculated when scale does not equal 'none'")
    }

    # check if individual specific allele frequencies are needed
    if(freq.type == "individual"){
        if(is.null(pcMat)){
            stop("pcMat is required to use individual-specific allele frequencies")
        }
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

    # check if population allele frequencies are needed
    if(freq.type == "population"){
        if(training.set$n == scan.include$n){
            if(verbose){ message(paste("Using all",training.set$n,"Samples to Estimate Population Allele Frequencies")) }
        }else{
            if(verbose){ message(paste("Using",training.set$n,"Samples in training.set to Estimate Population Allele Frequencies")) }
        }
    }


    # set up GDS files
    if(write.to.gds){
        # if no gds.prefix specified
        if(is.null(gds.prefix)){ gds.prefix <- "tmp" }

        if(verbose){ message("Creating GDS file for Allele Frequency Estimates:")}
        # create an empty GDS file
        gds <- createfn.gds(paste0(gds.prefix,"_freq.gds"))
        # add sample.id (closed)
        add.gdsn(gds, "sample.id", scan.include$value, compress="ZIP.max", closezip=TRUE)
        # add SNP info (open to append data)
        add.gdsn(gds, "snp.id", storage=class(snp.include$value), valdim=0, compress="ZIP.max")
        add.gdsn(gds, "snp.position", storage="integer", valdim=0, compress="ZIP.max")
        add.gdsn(gds, "snp.chromosome", storage="uint8", valdim=0, compress="ZIP.max")

        if(freq.type == "population"){
            add.gdsn(gds, "pop.freq", storage="float32", valdim=0, compress="ZIP.max")
        }
        if(freq.type == "individual"){
            # create node for individual-specific allele frequency data
            isaf.node <- add.gdsn(gds, "genotype", storage="float32", valdim=c(scan.include$n, 0), compress="ZIP.max")
            put.attr.gdsn(isaf.node, "sample.order")
        }
        sync.gds(gds)

        if(verbose){ message("Estimating Allele Frequencies...") }
        # loop through SNP blocks
        for(bn in 1:snp.blocks$n){
            # keep track of time for rate reporting
            startTime <- Sys.time()

            # index of SNPs for this block
            snp.block.idx <- snp.blocks$start[bn]:snp.blocks$end[bn]

            # append SNP info
            append.gdsn(index.gdsn(gds, "snp.id"), val = snp.include$value[snp.block.idx])
            if(class(genoData) == "GenotypeData"){
                append.gdsn(index.gdsn(gds, "snp.position"), val = getPosition(genoData, index = snp.include$index[snp.block.idx]))
                append.gdsn(index.gdsn(gds, "snp.chromosome"), val = getChromosome(genoData, index = snp.include$index[snp.block.idx]))
                # load genotype data
                geno <- getGenotypeSelection(genoData, snp = snp.include$index[snp.block.idx], scan = training.set$index, drop = FALSE)
            }else if(class(genoData) == "SeqVarData"){
                # set a filter to this SNP block
                seqSetFilter(genoData, variant.sel = snp.include$index[snp.block.idx], sample.sel = training.set$index, verbose = FALSE)
                append.gdsn(index.gdsn(gds, "snp.position"), val = seqGetData(genoData, "position"))
                append.gdsn(index.gdsn(gds, "snp.chromosome"), val = seqGetData(genoData, "chromosome"))
                # load genotype data
                geno <- t(altDosage(genoData))
            }

            if(freq.type == "population"){
                # estimate population allele frequencies
                phat <- 0.5*rowMeans(geno, na.rm=TRUE)
                append.gdsn(index.gdsn(gds, "pop.freq"), val = phat)
            }

            if(freq.type == "individual"){
                # estimate individual-specific allele frequencies
                muhat <- estISAF(geno, pcMat, pcProd, transpose = TRUE)
                append.gdsn(index.gdsn(gds, "genotype"), val=muhat)
            }
            sync.gds(gds)

            endTime <- Sys.time()
            rate <- format(endTime - startTime, digits=4)
            if(verbose){ message(paste("...SNP Block",bn,"of",snp.blocks$n,"Completed -", rate)) }
        }

        # put gds in readmode
        readmode.gdsn(index.gdsn(gds, "snp.id"))
        readmode.gdsn(index.gdsn(gds, "snp.position"))
        readmode.gdsn(index.gdsn(gds, "snp.chromosome"))
        if(freq.type == "population"){ readmode.gdsn(index.gdsn(gds, "pop.freq")) }
        if(freq.type == "individual"){ readmode.gdsn(isaf.node) }
        # cleanup
        filename <- gds$filename
        closefn.gds(gds)
        cleanup.gds(filename, verbose=FALSE)

        # read in created GDS
        freqData <- GenotypeData(GdsGenotypeReader(filename))

        # SNPs to include in the analysis (from freqData)
        snp.include.freq <- getSnpIndex(freqData, snp.include$value, chromosome)
        # Scans to include in the analysis (from freqData)
        scan.include.freq <- getScanIndex(freqData, scan.include$value)
        

    	if(verbose){ message("Creating GDS file for PC-Relate Results:")}    	
    	# create an empty GDS file for results
    	gds <- createfn.gds(paste0(gds.prefix,"_pcrelate.gds"))
    	# add sample.id (closed)
    	add.gdsn(gds, "sample.id", scan.include$value, compress="ZIP.max", closezip=TRUE)
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
                if(scan.blocks$n == 1){ message("Computing Relatedness Estimates...") }
                else{ message("Computing Relatedness Estimates for Scan Block Pair (",j,",",i,")...") }
            }
            # index of scans for block i
    		i.block.idx <- scan.blocks$start[i]:scan.blocks$end[i]
    		i.block.n <- length(i.block.idx)


    		# set up empty matrix for kinship coef
    		kin <- matrix(0, nrow=i.block.n, ncol=j.block.n)
            if(scale == "overall"){
                kinscale <- matrix(0, nrow=i.block.n, ncol=j.block.n)
            }
            # set up empty matrix for IBD probs
    		if(ibd.probs){
    			k2 <- matrix(0, nrow=i.block.n, ncol=j.block.n)
                k0 <- matrix(0, nrow=i.block.n, ncol=j.block.n)
                if(scale == "overall"){
                    if(freq.type == "individual"){
                        k2scale <- matrix(0, nrow=i.block.n, ncol=j.block.n)
                        k0scale <- matrix(0, nrow=i.block.n, ncol=j.block.n)
                    }
                    if(freq.type == "population"){
                        kscale <- matrix(0, nrow=i.block.n, ncol=j.block.n)
                    }
                }
    		}
            # set up empty matrix for number of SNPs used
    		nsnp <- matrix(0, nrow=i.block.n, ncol=j.block.n)


    		# loop through SNP blocks
    		for(bn in 1:snp.blocks$n){
    			# keep track of time for rate reporting
    			startTime <- Sys.time()

    			# index of SNPs for this block
    			snp.block.idx <- snp.blocks$start[bn]:snp.blocks$end[bn]

                # load genotype data
                if(class(genoData) == "GenotypeData"){
                    geno <- getGenotypeSelection(genoData, snp = snp.include$index[snp.block.idx], scan = scan.include$index[unique(c(i.block.idx,j.block.idx))], drop = FALSE)
                }else if(class(genoData) == "SeqVarData"){
                    # set a filter to these scans and this variant block
                    seqSetFilter(genoData, variant.sel = snp.include$index[snp.block.idx], sample.sel = scan.include$index[unique(c(i.block.idx,j.block.idx))], verbose = FALSE)
                    geno <- t(altDosage(genoData))
                }                
                # index for which columns of geno are in each block                 			
    			i.geno.idx <- which(colnames(geno) %in% scan.include$value[i.block.idx])
    			j.geno.idx <- which(colnames(geno) %in% scan.include$value[j.block.idx])

                # get population allele frequencies
                if(freq.type == "population"){
                    if(write.to.gds){
                        # load saved frequencies
                        phat <- getVariable(freqData, "pop.freq", index = snp.include.freq$index[snp.block.idx])
                    }
                    else{
                        # estimate frequencies
                        phat <- 0.5*rowMeans(geno, na.rm=TRUE)
                    }                    

                    # which SNPs to filter (pop freq missing or too big/small)
                    pop.filt <- (is.na(phat) | phat == 0 | phat == 1 | phat < MAF | phat > (1-MAF))
                    phat <- phat[!pop.filt]
                    geno <- geno[!pop.filt,]
                    rm(pop.filt)

                    # opposite allele frequency
                    if(scale != "none"){
                        qhat <- 1 - phat
                        phatqhat <- phat*qhat
                    }
                }

                # get individual-specific allele frequencies
                if(freq.type == "individual"){                    
                    if(write.to.gds){
                        # load saved individual specific allele frequency data
                        muhat <- getGenotypeSelection(freqData, snp = snp.include.freq$index[snp.block.idx], scan = scan.include.freq$index[unique(c(i.block.idx,j.block.idx))], drop = FALSE)
                    }else{
                        # estimate individual-specific allele frequencies
                        muhat <- estISAF(geno[,scan.include$value %in% training.set$value], pcMat, pcProd, transpose = FALSE)
                    }

                    # which SNPs to filter (isaf missing or too big/small)
                    isaf.filt <- which(is.na(muhat) | muhat == 0 | muhat == 1 | muhat < MAF | muhat > (1-MAF))
                    muhat[isaf.filt] <- NA
                    geno[isaf.filt] <- NA
                    rm(isaf.filt)

                    # opposite allele frequency
                    if(scale != "none"){
                        if(ibd.probs){
                            quhat <- 1-muhat
                            muhatquhat <- muhat*quhat
                        }else{
                            muhatquhat <- muhat*(1-muhat)
                        }
                    }
                }

                # determine the number of snps used for each pair
                nonmiss <- !is.na(geno)
                nsnp <- nsnp + crossprod(nonmiss[,i.geno.idx], nonmiss[,j.geno.idx])
                # index of SNPs to filter
                filt.idx <- which(!nonmiss); rm(nonmiss)


                # kinship
                if(freq.type == "population"){
                    Z <- geno - 2*phat
                    if(scale == "variant"){
                        sigma <- sqrt(phatqhat)
                        Z <- Z/sigma
                    }
                }
                if(freq.type == "individual"){
                    Z <- geno - 2*muhat
                    if(scale == "variant"){
                        Z <- Z/sqrt(muhatquhat)
                    }
                }
                # set missing values to 0 (so no contribution)
                Z[filt.idx] <- 0
                # update kin
                kin <- kin + crossprod(Z[,i.geno.idx], Z[,j.geno.idx]); rm(Z)

                # overall scaling             
                if(scale == "overall"){
                    if(freq.type == "population"){
                        phatqhatMat <- matrix(sqrt(phatqhat), nrow=length(phatqhat), ncol=ncol(geno))
                        # set missing values to 0 (so no contribution)
                        phatqhatMat[filt.idx] <- 0
                        kinscale <- kinscale + crossprod(phatqhatMat[,i.geno.idx], phatqhatMat[,j.geno.idx])
                        rm(phatqhatMat)                        
                    }
                    if(freq.type == "individual"){
                        # set missing values to 0 (so no contribution)
                        muhatquhat[filt.idx] <- 0
                        kinscale <- kinscale + crossprod(sqrt(muhatquhat[,i.geno.idx]), sqrt(muhatquhat[,j.geno.idx]))
                        if(!ibd.probs){ rm(muhatquhat) }
                    }
                }


                if(ibd.probs){
                    # k2
                    # dominance coded matrix 
                    Zd <- matrix(0, nrow = nrow(geno), ncol = ncol(geno))
                    if(freq.type == "population"){
                        idx0 <- which(geno == 0)
                        idx0p <- idx0 %% nrow(Zd); idx0p[idx0p == 0] <- nrow(Zd)
                        Zd[idx0] <- phat[idx0p]; rm(idx0); rm(idx0p)
                        idx2 <- which(geno == 2)
                        idx2p <- idx2 %% nrow(Zd); idx2p[idx2p == 0] <- nrow(Zd)
                        Zd[idx2] <- qhat[idx2p]; rm(idx2); rm(idx2p)
                        Zd <- Zd - phatqhat
                        if(scale == "variant"){
                            Zd <- Zd/phatqhat
                        }
                    }                    
                    if(freq.type == "individual"){
                        idx0 <- which(geno == 0)
                        Zd[idx0] <- muhat[idx0]; rm(idx0)
                        idx2 <- which(geno == 2)
                        Zd[idx2] <- quhat[idx2]; rm(idx2)                    
                        Zd <- Zd - muhatquhat
                        if(scale == "variant"){
                            Zd <- Zd/muhatquhat
                        }
                    }
                    # set missing values to 0 (so no contribution)
                    Zd[filt.idx] <- 0
                    # update k2
                    k2 <- k2 + crossprod(Zd[,i.geno.idx], Zd[,j.geno.idx]); rm(Zd)

                    # k0
                    # indicator matrices of homozygotes
                    IAA <- geno == 2
                    Iaa <- geno == 0
                    if(scale == "variant"){
                        if(freq.type == "population"){                        
                            IAA <- IAA/phat^2
                            Iaa <- 0.5*Iaa/qhat^2
                            # factor of 1/2 for IAA*Iaa comes in from splitting up the two terms AA,aa vs aa,AA
                        }
                        if(freq.type == "individual"){
                            IAA <- IAA/muhat^2
                            Iaa <- 0.5*Iaa/quhat^2                            
                            rm(quhat)
                        }
                    }
                    # set missing values to 0 (so no contribution)
                    IAA[filt.idx] <- 0
                    Iaa[filt.idx] <- 0
                    # update k0
                    k0 <- k0 + crossprod(IAA[,i.geno.idx], Iaa[,j.geno.idx]) + crossprod(Iaa[,i.geno.idx], IAA[,j.geno.idx]); rm(IAA); rm(Iaa)

                    # overall scaling
                    if(scale == "overall"){
                        if(freq.type == "population"){
                            phatqhatMat <- matrix(phatqhat, nrow=length(phatqhat), ncol=ncol(geno))
                            # set missing values to 0 (so no contribution)
                            phatqhatMat[filt.idx] <- 0
                            kscale <- kscale + crossprod(phatqhatMat[,i.geno.idx], phatqhatMat[,j.geno.idx]); rm(phatqhatMat)
                            # this is k2scale; 2*kscale is k0scale
                        }
                        if(freq.type == "individual"){                            
                            k2scale <- k2scale + crossprod(muhatquhat[,i.geno.idx], muhatquhat[,j.geno.idx])
                            # set missing values to 0 (so no contribution)
                            muhat[filt.idx] <- 0
                            quhat[filt.idx] <- 0
                            k0scale <- k0scale + crossprod(muhat[,i.geno.idx]^2, quhat[,j.geno.idx]^2) + crossprod(quhat[,i.geno.idx]^2, muhat[,j.geno.idx]^2)
                        }                        
                    }
                }

    			endTime <- Sys.time()
    			rate <- format(endTime - startTime, digits=4)
    			if(verbose){ message(paste("...SNP Block",bn,"of",snp.blocks$n,"Completed -", rate)) }
    		}


    		# compute kinship/inbreeding estimates
            if(scale == "overall"){
                kin <- kin/(4*kinscale); rm(kinscale)                
            }else{
                kin <- kin/(4*nsnp)
            }


            # small sample correction
            kincorrect <- NULL
            if(scan.blocks$n == 1 & correct & freq.type == "individual" & scale != "none"){
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
                diag(kin) <- 0.5*(1 + f.tmp); rm(f.tmp)
            }

            k2correct <- NULL
            if(ibd.probs){
                # compute k2 estimates
                if(scale == "overall"){
                    if(freq.type == "population"){
                        k2 <- k2/kscale
                    }
                    if(freq.type == "individual"){
                        k2 <- k2/k2scale; rm(k2scale)
                    }                    
                }else{
                    k2 <- k2/nsnp
                }
                # correction for HW departure                
                if(i == j){ 
                    f.j <- 2*diag(kin) - 1
                    fprod <- tcrossprod(f.j)
                }else{
                    f.i <- 2*diag(read.gdsn(kinship.node, start=c(i.block.idx[1], i.block.idx[1]), count=c(i.block.n, i.block.n))) - 1
                    fprod <- tcrossprod(f.i,f.j)
                }
                k2 <- k2 - fprod

                # small sample correction
                if(scan.blocks$n == 1 & correct & freq.type == "individual" & scale != "none"){                
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
                if(scale == "overall"){
                    if(freq.type == "population"){
                        k0 <- k0/(2*kscale); rm(kscale)
                    }
                    if(freq.type == "individual"){
                        k0 <- k0/k0scale; rm(k0scale)
                    }
                }else{
                    k0 <- k0/nsnp
                }

                # index for 1st deg rels
                not1stDegidx <- which(kin < 2^(-5/2))
                k0[not1stDegidx] <- 1 - 4*kin[not1stDegidx] + k2[not1stDegidx]               
            }


            # set up output
            if(i == j){
                nsnp[upper.tri(nsnp)] <- NA
                if(ibd.probs){
                    k2[lower.tri(k2)] <- k0[lower.tri(k0)]                    
                    diag(k2) <- NA
                }
            }

            # write output
            if(write.to.gds){
                write.gdsn(kinship.node, kin, start=c(i.block.idx[1], j.block.idx[1]), count=c(i.block.n, j.block.n))
                write.gdsn(nsnp.node, t(nsnp), start=c(j.block.idx[1], i.block.idx[1]), count=c(j.block.n, i.block.n))
                if(ibd.probs){
                    write.gdsn(ibd.node, k2, start=c(i.block.idx[1], j.block.idx[1]), count=c(i.block.n, j.block.n))
                }
                if(i != j){
                    write.gdsn(kinship.node, t(kin), start=c(j.block.idx[1], i.block.idx[1]), count=c(j.block.n, i.block.n))
                    write.gdsn(nsnp.node, matrix(NA, nrow=i.block.n, ncol=j.block.n), start=c(i.block.idx[1], j.block.idx[1]), count=c(i.block.n, j.block.n))                 
                    if(ibd.probs){
                        write.gdsn(ibd.node, t(k0), start=c(j.block.idx[1], i.block.idx[1]), count=c(j.block.n, i.block.n))
                    }
                }
                sync.gds(gds)
            }       

    	} # j loop
    } # i loop

    if(write.to.gds){
        # add additional information        
        add.gdsn(gds, "kincorrect", kincorrect, closezip = TRUE)
        add.gdsn(gds, "k2correct", k2correct, closezip = TRUE)
        add.gdsn(gds, "freq.type", freq.type, storage="character", closezip=TRUE)
        add.gdsn(gds, "scale", scale, storage="character", closezip=TRUE)        
    	# compress nodes and put in readmode
    	compression.gdsn(kinship.node, compress="ZIP.max")    	
    	compression.gdsn(nsnp.node, compress="ZIP.max")
    	readmode.gdsn(kinship.node)    	
    	readmode.gdsn(nsnp.node)
    	if(ibd.probs){
    		compression.gdsn(ibd.node, compress="ZIP.max")
    		readmode.gdsn(ibd.node)
    	}
    	# cleanup
    	filename <- gds$filename
    	closefn.gds(gds)
    	close(freqData)
    	cleanup.gds(filename, verbose=FALSE)

        out <- paste("Output saved to", filename, sep=" ")

    }else{
    	if(ibd.probs){
    		out <- list(sample.id = scan.include$value, kinship = kin, ibd.probs = k2, nsnp = nsnp, kincorrect = kincorrect, k2correct = k2correct, call = match.call(), freq.type = freq.type, scale = scale)
    	}else{
    		out <- list(sample.id = scan.include$value, kinship = kin, nsnp = nsnp, kincorrect = kincorrect, call = match.call(), freq.type = freq.type, scale = scale)
    	}
      class(out) <- "pcrelate"
    }

    if(class(genoData) == "SeqVarData"){ seqSetFilter(genoData, sample.sel = seqFilt.original$sample.sel, variant.sel = seqFilt.original$variant.sel, verbose = FALSE) }

    return(out)
}



