pcrelate2 <- function(	gdsobj,
						pcs,
						scale = c('overall', 'variant', 'none'),
						ibd.probs = TRUE,
						sample.include = NULL,
						training.set = NULL,
						sample.block.size = 5000,
						## snp.include = NULL,
						## snp.block.size = 10000,
						maf.thresh = 0.01,
						maf.bound.method = c('filter', 'truncate'),
						small.samp.correct = FALSE,
						num.cores = 1,
						verbose = TRUE){

	# checks
        scale <- match.arg(scale)
        maf.bound.method <- match.arg(maf.bound.method)
        if(is.null(sample.include)) sample.include <- .readSampleId(gdsobj)
        .pcrelateChecks(pcs = pcs, scale = scale, ibd.probs = ibd.probs, sample.include = sample.include, training.set = training.set, 
					maf.thresh = maf.thresh)
	
	# set up number of cores
	sys.cores <- parallel::detectCores(logical = TRUE)
	doMC::registerDoMC(cores = min(c(num.cores, sys.cores)))
	if(verbose) message('Using ', min(c(num.cores, sys.cores)), ' CPU cores')

        # number of sample blocks
        nsampblock <- ceiling(length(sample.include)/sample.block.size)
	# list of samples in each block
	if(nsampblock > 1){
		samp.blocks <- unname(split(sample.include, cut(1:length(sample.include), nsampblock)))
	}else{
		samp.blocks <- list(sample.include)
	}

	if(nsampblock == 1){
		if(verbose) message(length(sample.include), ' samples to be included in the analysis...')

		# create matrix of PCs
		V <- .createPCMatrix(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include)

		# matrix product of V
		VVtVi <- .calcISAFBetaPCProd(V = V, training.set = training.set, verbose = verbose)

		### Stephanie to build in variant iterators here ###
		# number of snp blocks
		## nsnpblock <- ceiling(length(snp.include)/snp.block.size)
		## # list of snps in each block
		## if(nsnpblock > 1){
		## 	snp.blocks <- unname(split(snp.include, cut(1:length(snp.include), nsnpblock)))
		## }else{
		## 	snp.blocks <- list(snp.include)
		## }
                snp.blocks <- .snpBlocks(gdsobj)
                nsnpblock <- length(snp.blocks)
		if(verbose) message('Running PC-Relate analysis for ', length(sample.include), ' samples using ', length(unlist(snp.blocks)), ' SNPs in ', nsnpblock, ' blocks...')

		# for each snp block
		matList <- foreach(k = 1:nsnpblock, .combine = .matListCombine, .inorder = FALSE, .multicombine = FALSE) %dopar% {
			# read genotype data for the block
			#seqSetFilter(gdsobj, variant.id = snp.blocks[[k]])
                        #G <- altDosage(gdsobj)
                        G <- .readGeno(gdsobj, sample.include, snp.index = snp.blocks[[k]])

			# calculate ISAF betas
			beta <- .calcISAFBeta(G = G, VVtVi = VVtVi)

			# calculate PC-Relate estimates
			matList.block <- .pcrelateVarBlock(G = G, beta = beta, V = V, idx = 1:nrow(V), jdx = 1:nrow(V), scale = scale, ibd.probs = ibd.probs, maf.thresh = maf.thresh, maf.bound.method = maf.bound.method)
			matList.block
		}

		# take ratios to compute final estimates
		estList <- .pcrelateCalcRatio(matList = matList, scale = scale, ibd.probs = ibd.probs)

		# cast to data.tables
		kinSelf <- data.table(ID = rownames(estList$kin), f = 2*diag(estList$kin) - 1, nsnp = diag(estList$nsnp))
		kinBtwn <- .estListToDT(estList, drop.lower = TRUE)


	}else if(nsampblock > 1){
		if(verbose) message(length(sample.include), ' samples to be included in the analysis, split into ', nsampblock, ' blocks...')

		# calculate betas for individual specific allele frequencies
		beta <- calcISAFBeta(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include, training.set = training.set, verbose = verbose)		
		### beta is a matrix of variants x PCs; needs to be saved, not sure of the best format or the best way to split this up ###
		
		# compute estimates for current (pair of) sample block(s)
		### this is where we would parallelize with multiple jobs ###
		kinSelf <- NULL
		kinBtwn <- NULL
		for(i in 1:nsampblock){
			for(j in i:nsampblock){
				if(verbose) message('Running PC-Relate analysis for sample block pair (', i, ',', j, ')')
				# compute estimates for the (pair of) sample block(s)
				tmp <- pcrelateSampBlock(	gdsobj = gdsobj, pcs = pcs, betaobj = beta, sample.include.block1 = samp.blocks[[i]], sample.include.block2 = samp.blocks[[j]],
											scale = scale, ibd.probs = ibd.probs,
											maf.thresh = maf.thresh, maf.bound.method = maf.bound.method, verbose = verbose)

				# update results with this (pair of) sample block(s)
				if(i == j) kinSelf <- rbind(kinSelf, tmp$kinSelf)
				kinBtwn <- rbind(kinBtwn, tmp$kinBtwn)
			}
		}
	}

	# order samples
	setkey(kinSelf, ID)
	setkey(kinBtwn, ID1, ID2)

	### post processing after putting all samples together
	if(ibd.probs){
		# correct k2 for HW departure and use alternate k0 estimator for non-1st degree relatives
		kinBtwn <- .fixIBDEst(kinBtwn = kinBtwn, kinSelf = kinSelf)
	}

	### need to rewrite the small sample correction at some point - take care of this later ###


	return(list(kinBtwn = kinBtwn, kinSelf = kinSelf))
}





.pcrelateChecks <- function(pcs, scale, ibd.probs, sample.include, training.set, maf.thresh){
	# check parameters
	if(scale == 'none' & ibd.probs) stop('`ibd.probs` must be FALSE when `scale` == none')
	if(maf.thresh < 0 | maf.thresh > 0.5) stop('maf.thresh must be in [0,0.5]')
	# check training.set
	if(!is.null(training.set) & !all(training.set %in% sample.include)) stop('All samples in training.set must be in sample.include')
	# check pcs
	if(class(pcs) != 'matrix' | is.null(rownames(pcs))) stop('pcs should be a matrix of PCs with rownames set to sample.ids')
	if(!all(sample.include %in% rownames(pcs))) stop('All samples in sample.include must be in pcs')
}


### function to match samples and create PC matrix
.createPCMatrix <- function(gdsobj, pcs, sample.include){
    # set filter to sample.include
    #seqSetFilter(gdsobj, sample.id = sample.include)
    sample.id <- intersect(.readSampleId(gdsobj), sample.include)

    # subset and re-order pcs if needed
    V <- pcs[match(sample.id, rownames(pcs)), , drop = FALSE]
    # append intercept
    V <- cbind(rep(1, nrow(V)), V)
    return(V)
}

# function to get 
.calcISAFBetaPCProd <- function(V, training.set, verbose = TRUE){
	if(!is.null(training.set)){
		idx <- rownames(V) %in% training.set
		VVtVi <- tcrossprod(V[idx,], chol2inv(chol(crossprod(V[idx,]))))
		if(verbose) message('Betas for ', ncol(V) - 1, ' PC(s) will be calculated using ', sum(idx), ' samples in training.set...')
	}else{
		idx <- NULL
		VVtVi <- tcrossprod(V, chol2inv(chol(crossprod(V))))
		if(verbose) message('Betas for ', ncol(V) - 1, ' PC(s) will be calculated using all ', nrow(V), ' samples in sample.include...')
	}
	return(list(val = VVtVi, idx = idx))
}

# function to do actual calculation of betas
.calcISAFBeta <- function(G, VVtVi){
	# impute missing genotype values
        G <- .meanImpute(G, freq=0.5*colMeans(G, na.rm=TRUE))

	# calculate beta
	if(is.null(VVtVi$idx)){
		beta <- crossprod(G, VVtVi$val)
	}else{
		beta <- crossprod(G[VVtVi$idx, ], VVtVi$val)
	}
	return(beta)
}



### this function does the pcrelate estimation for a variant block
.pcrelateVarBlock <- function(G, beta, V, idx, jdx, scale, ibd.probs, maf.thresh, maf.bound.method){

	# make sure order of G, beta, and V all line up
	if(!all.equal(colnames(G), rownames(beta))) stop('G and beta do not match')
	if(!all.equal(rownames(G), rownames(V))) stop('G and V do not match')

	# estimate individual specific allele frequencies
	mu <- .estISAF(beta = beta, V = V, bound.thresh = maf.thresh, bound.method = maf.bound.method)

	# compute values for estimates
	if(scale == 'overall'){
		matList <- .pcrCalcOvr(G = G, mu = mu, idx = idx, jdx = jdx, ibd.probs = ibd.probs)
	}else if(scale == 'variant'){
		matList <- .pcrCalcVar(G = G, mu = mu, idx = idx, jdx = jdx, ibd.probs = ibd.probs)
	}else if(scale == 'none'){
		matList <- .pcrCalcNone(G = G, mu = mu, idx = idx, jdx = jdx)
	}

	return(matList)
}

.matListCombine <- function(...){
    mapply(FUN = "+", ..., SIMPLIFY = FALSE)
}


### these functions are for estimating individual specific allele frequencies
.estISAF <- function(beta, V, bound.thresh, bound.method){
	# get ISAF estimates (i.e. 0.5*fitted values)
	mu <- 0.5*tcrossprod(V, beta)
	# fix boundary cases
	if(bound.method == 'truncate'){
		mu <- apply(mu, 2, .muBoundaryTrunc, thresh = bound.thresh)
	}else if(bound.method == 'filter'){
		mu <- apply(mu, 2, .muBoundaryFilt, thresh = bound.thresh)
	}else{
		stop('bound.method must be one of `truncate` or `filter`')
	}
	return(mu)
}

# function to truncate boundary issues with ISAFs
.muBoundaryTrunc <- function(x, thresh){
    x[x < thresh] <- thresh
    x[x > (1 - thresh)] <- (1 - thresh)
    return(x)
}

# function to filter boundary issues with ISAFs
.muBoundaryFilt <- function(x, thresh){
    x[x < thresh] <- NA
    x[x > (1 - thresh)] <- NA
    return(x)
}



### all of the functions below are used for computing kinship and ibd estimates ###
.pcrCalcOvr <- function(G, mu, ibd.probs, idx, jdx){
	# compute mu(1 - mu)
	muqu <- mu*(1 - mu)

	# compute number of observed snps by pair
	nonmiss <- !(is.na(G) | is.na(mu))
	nsnp <- tcrossprod(nonmiss)

	# index of missing values
	filt.idx <- which(!nonmiss)

	# compute kinship values
	kinList <- .pcrCalcKinOvr(G, mu, muqu, filt.idx, idx, jdx)

	if(ibd.probs){
		# compute ibd values
		ibdList <- .pcrCalcIBDOvr(G, mu, muqu, nonmiss, filt.idx, idx, jdx)

		return(list(kinNum = kinList$kinNum, 
					kinDen = kinList$kinDen,
					k0Num = ibdList$k0Num,
					k0Den = ibdList$k0Den,
					k2Num = ibdList$k2Num,
					k2Den = ibdList$k2Den,
					nsnp = nsnp))

	}else{
		return(list(kinNum = kinList$kinNum,
					kinDen = kinList$kinDen,
					nsnp = nsnp))
	}
}

.pcrCalcKinOvr <- function(G, mu, muqu, filt.idx, idx, jdx){
	# residuals
	R <- G - 2*mu

	# set missing values to 0 (i.e. no contribution from that variant)
	R[filt.idx] <- 0
	muqu[filt.idx] <- 0

	# numerator (crossprod of residuals)
	kinNum <- tcrossprod(R[idx,], R[jdx,])
	# denominator
	kinDen <- tcrossprod(sqrt(muqu[idx,]), sqrt(muqu[jdx,]))

	return(list(kinNum = kinNum, kinDen = kinDen))
}

.pcrCalcIBDOvr <- function(G, mu, muqu, nonmiss, filt.idx, idx, jdx){
	# opposite allele frequency
	qu <- 1 - mu

    # indicator matrices of homozygotes
	Iaa <- G == 0 & nonmiss
	IAA <- G == 2 & nonmiss

	# dominance coded genotype matrix
	Gd <- mu
	Gd[G == 1 & nonmiss] <- 0
	Gd[IAA] <- qu[IAA]
	Gd <- Gd - muqu

	# set missing values to 0 (i.e. no contribution from that variant)
	Iaa[filt.idx] <- 0
	IAA[filt.idx] <- 0
	Gd[filt.idx] <- 0
	mu[filt.idx] <- 0
	qu[filt.idx] <- 0
	muqu[filt.idx] <- 0

	# k0	
	k0Num <- tcrossprod(IAA[idx,], Iaa[jdx,]) + tcrossprod(Iaa[idx,], IAA[jdx,])
	k0Den <- tcrossprod(mu[idx,]^2, qu[jdx,]^2) + tcrossprod(qu[idx,]^2, mu[jdx,]^2)

	# k2
	k2Num <- tcrossprod(Gd[idx,], Gd[jdx,])
	k2Den <- tcrossprod(muqu[idx,], muqu[jdx,])

	return(list(k0Num = k0Num, k0Den = k0Den, k2Num = k2Num, k2Den = k2Den))
}


.pcrCalcVar <- function(G, mu, ibd.probs, idx, jdx){
	# compute mu(1 - mu)
	muqu <- mu*(1 - mu)

	# compute number of observed snps by pair
	nonmiss <- !(is.na(G) | is.na(mu))
	nsnp <- tcrossprod(nonmiss)

	# index of missing values
	filt.idx <- which(!nonmiss)

	# compute kinship values
	kinNum <- .pcrCalcKinVar(G, mu, muqu, filt.idx, idx, jdx)

	if(ibd.probs){
		# compute ibd values
		ibdList <- .pcrCalcIBDVar(G, mu, muqu, nonmiss, filt.idx, idx, jdx)

		return(list(kinNum = kinNum,
					k0Num = ibdList$k0Num,
					k2Num = ibdList$k2Num,
					nsnp = nsnp))

	}else{
		return(list(kinNum = kinNum, nsnp = nsnp))
	}
}

.pcrCalcKinVar <- function(G, mu, muqu, filt.idx, idx, jdx){
	# scaled residuals
	R <- (G - 2*mu)/sqrt(muqu)

	# set missing values to 0 (i.e. no contribution from that variant)
	R[filt.idx] <- 0

	# numerator (crossprod of scaled residuals)
	kinNum <- tcrossprod(R[idx,], R[jdx,])
	
	return(kinNum)
}

.pcrCalcIBDVar <- function(G, mu, muqu, nonmiss, filt.idx, idx, jdx){
	# opposite allele frequency
	qu <- 1 - mu

	# indicator matrices of homozygotes
	Iaa <- G == 0 & nonmiss
	IAA <- G == 2 & nonmiss

	# dominance coded genotype matrix
	Gd <- mu
	Gd[G == 1 & nonmiss] <- 0
	Gd[IAA] <- qu[IAA] 

	# scale 
	Iaa <- 0.5*Iaa/qu^2 # factor of 1/2 for IAA*Iaa comes in from splitting up the two terms AA,aa vs aa,AA
	IAA <- IAA/mu^2
	Gd <- Gd/muqu - 1

	# set missing values to 0 (i.e. no contribution from that variant)
	Iaa[filt.idx] <- 0
	IAA[filt.idx] <- 0
	Gd[filt.idx] <- 0

	# k0
	k0Num <- tcrossprod(IAA[idx,], Iaa[jdx,]) + tcrossprod(Iaa[idx,], IAA[jdx,])

	# k2
	k2Num <- tcrossprod(Gd[idx,], Gd[jdx,])

	return(list(k0Num = k0Num, k2Num = k2Num))
}


.pcrCalcNone <- function(G, mu, idx, jdx){
	# compute number of observed snps by pair
	nonmiss <- !(is.na(G) | is.na(mu))
	nsnp <- tcrossprod(nonmiss)

	# index of missing values
	filt.idx <- which(!nonmiss)
	rm(nonmiss)

	# compute kinship values
	kin <- .pcrCalcKinNone(G, mu, filt.idx, idx, jdx)
	return(list(kin = kin, nsnp = nsnp))
}

.pcrCalcKinNone <- function(G, mu, filt.idx, idx, jdx){
	# residuals
	R <- G - 2*mu

	# set missing values to 0 (i.e. no contribution from that variant)
	R[filt.idx] <- 0

	# crossprod of scaled residuals
	kin <- tcrossprod(R[idx,], R[jdx,])

	return(kin)
}


### functions for post processing on variant block level
.pcrelateCalcRatio <- function(matList, scale, ibd.probs){
	# compute final estimates
	if(scale == 'overall'){
		kin <- matList$kinNum/(4*matList$kinDen)
	}else if(scale == 'variant'){
		kin <- matList$kin/(4*matList$nsnp)
	}

	if(ibd.probs){
		if(scale == 'overall'){
			k2 <- matList$k2Num/matList$k2Den
			k0 <- matList$k0Num/matList$k0Den
		}else if(scale == 'variant'){
			k2 <- matList$k2/matList$nsnp
			k0 <- matList$k0/matList$nsnp
		}
		return(list(kin = kin, k0 = k0, k2 = k2, nsnp = matList$nsnp))

	}else{
		return(list(kin = kin, nsnp = matList$nsnp))
	}
}

		
### functions for post processing on sample block level
.estListToDT <- function(x, drop.lower){
	estDT <- lapply(x, .estMelt, drop.lower = drop.lower)
	for(k in 1:length(estDT)){
		setnames(estDT[[k]], 'value', names(estDT)[k])
	}
	# merge those data.tables into one data.table
	estDT <- Reduce(merge, estDT)
	return(estDT)
}

.estMelt <- function(x, drop.lower){
	if(drop.lower){
        x[lower.tri(x, diag = TRUE)] <- NA
    }
    x <- as.data.table(reshape2::melt(x, varnames = c('ID1', 'ID2'), na.rm = TRUE))
    setkey(x, ID1, ID2)
}


### functions for final processing
.fixIBDEst <- function(kinBtwn, kinSelf){
    # correct k2 for HW departure
    kinBtwn <- merge(kinBtwn, kinSelf[,.(ID, f)], by.x = 'ID2', by.y = 'ID')
    setnames(kinBtwn, 'f', 'f.2')
    kinBtwn <- merge(kinBtwn, kinSelf[,.(ID, f)], by.x = 'ID1', by.y = 'ID')
    setnames(kinBtwn, 'f', 'f.1')
    kinBtwn[, k2 := k2 - f.1*f.2][, `:=`(f.1 = NULL, f.2 = NULL)]
    
    # use alternate k0 estimator for non-1st degree relatives
    kinBtwn[kin < 2^(-5/2), k0 := 1 - 4*kin + k2]
    
    return(kinBtwn)
}



### exported function for computing PC betas for individual specific allele frequency calculations ###
calcISAFBeta <- function(gdsobj, pcs, sample.include, training.set = NULL, snp.include, snp.block.size, verbose = TRUE){
	# checks - add some

	# create matrix of PCs
	V <- .createPCMatrix(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include)

	# matrix product of V
	VVtVi <- .calcISAFBetaPCProd(V = V, training.set = training.set, verbose = verbose)
	
	### Stephanie to build in variant iterators here ###
	# number of snp blocks
	## nsnpblock <- ceiling(length(snp.include)/snp.block.size)
	## # list of snps in each block
	## if(nsnpblock > 1){
	## 	snp.blocks <- unname(split(snp.include, cut(1:length(snp.include), nsnpblock)))
	## }else{
	## 	snp.blocks <- list(snp.include)
	## }
        snp.blocks <- .snpBlocks(gdsobj)
        nsnpblock <- length(snp.blocks)
	if(verbose) message('Calculating Indivdiual-Specific Allele Frequency betas for ', length(unlist(snp.blocks)), ' SNPs in ', nsnpblock, ' blocks...')

	beta <- foreach(k = 1:nsnpblock, .combine = rbind, .inorder = FALSE, .multicombine = TRUE) %dopar% {
		# read genotype data for the block
		#seqSetFilter(gdsobj, variant.id = snp.blocks[[k]])
		#G <- altDosage(gdsobj)
                G <- .readGeno(gdsobj, sample.include, snp.index = snp.blocks[[k]])

		# calculate ISAF betas
		beta.block <- .calcISAFBeta(G = G, VVtVi = VVtVi)
		beta.block
	}
	### rather than returning and rbinding here, we may want to be writing the output to something

	# non-parallel version
	# beta <- NULL
	# for(k in 1:nsnpblock){
	# 	# calculate betas for PCs at each variant
	# 	beta.block <- .calcISAFBeta(gdsobj = gdsobj, VVtVi = VVtVi, snp.include = snp.blocks[[k]])
	# 	if(verbose) message('...SNP Block ', k, ' of ', nsnpblock, ' Completed: ', nrow(beta.block), ' SNPs')
	# 	# rbind
	# 	beta <- rbind(beta, beta.block)
	# }

	return(beta)
}


### exported function that does the pcrelate estimation for a (pair of) sample block(s)
pcrelateSampBlock <- function(gdsobj, betaobj, pcs, sample.include.block1, sample.include.block2, scale, ibd.probs, maf.thresh, maf.bound.method, verbose = TRUE){

        # create (joint) PC matrix and indices
        sample.include <- unique(c(sample.include.block1, sample.include.block2))
	V <- .createPCMatrix(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include)
	idx <- which(rownames(V) %in% sample.include.block1)
	jdx <- which(rownames(V) %in% sample.include.block2)
	oneblock <- ifelse(all(idx == jdx), TRUE, FALSE)
	### slight inefficiency above because we create V for samples in block1 for each block2 when we don't have to if we are running serially; 
	### but this seems more straightforward to parallelize

	### Stephanie to build in variant iterators here ###
	# number of snp blocks
	## nsnpblock <- ceiling(length(snp.include)/snp.block.size)
	## # list of snps in each block
	## if(nsnpblock > 1){
	## 	snp.blocks <- unname(split(snp.include, cut(1:length(snp.include), nsnpblock)))
	## }else{
	## 	snp.blocks <- list(snp.include)
	## }
        snp.blocks <- .snpBlocks(gdsobj)
        nsnpblock <- length(snp.blocks)

	if(verbose) message('Running PC-Relate analysis using ', length(unlist(snp.blocks)), ' SNPs in ', nsnpblock, ' blocks...')
	# compute estimates for each variant block; sum along the way
	matList <- foreach(k = 1:nsnpblock, .combine = .matListCombine, .inorder = FALSE, .multicombine = FALSE) %dopar% {
		# read genotype data for the block
		#seqSetFilter(gdsobj, variant.id = snp.blocks[[k]])
		#G <- altDosage(gdsobj)
                G <- .readGeno(gdsobj, sample.include, snp.index = snp.blocks[[k]])

		# load betas for the current block of variants
		beta.block <- betaobj[rownames(betaobj) %in% snp.blocks[[k]], , drop = FALSE]
		### this line of code will probably be different if we save the betas; need to load correct betas

		# calculate PC-Relate estimates
		matList.block <- .pcrelateVarBlock(	G = G, beta = beta.block, V = V, idx = idx, jdx = jdx, scale = scale, ibd.probs = ibd.probs, maf.thresh = maf.thresh, maf.bound.method = maf.bound.method)
		matList.block
	}

	# non-parallel version
	# for(k in 1:nsnpblock){
	# 	# current variant block
	# 	tmpList <- .pcrelateVarBlock(	gdsobj = gdsobj, beta = beta, V = V, idx = idx, jdx = jdx, scale = scale, ibd.probs = ibd.probs, 
	# 									snp.include = snp.blocks[[k]], maf.thresh = maf.thresh, maf.bound.method = maf.bound.method)

	# 	# update overall (sum across variant blocks)
	# 	if(k == 1){
	# 		matList <- tmpList
	# 	}else{
	# 		for(m in 1:length(matList)){
	# 			matList[[m]] <- matList[[m]] + tmpList[[m]]
	# 		}
	# 	}
	# 	rm(tmpList)
	# 	if(verbose) message('...SNP Block ', k, ' of ', nsnpblock, ' Completed: ', length(snp.blocks[[k]]), ' SNPs')
	# }

	# take ratios to compute final estimates
	estList <- .pcrelateCalcRatio(matList = matList, scale = scale, ibd.probs = ibd.probs)

	# cast to data.tables
	if(oneblock){
		# self table
		kinSelf <- data.table(ID = rownames(estList$kin), f = 2*diag(estList$kin) - 1, nsnp = diag(estList$nsnp))
		setkey(kinSelf, ID)
		# between samples table
		kinBtwn <- .estListToDT(estList, drop.lower = TRUE)
	}else{
		kinSelf <- NULL
		# between samples table
		kinBtwn <- .estListToDT(estList, drop.lower = FALSE)
	}
	setkey(kinBtwn, ID1, ID2)
	
	return(list(kinSelf = kinSelf, kinBtwn = kinBtwn))
}
