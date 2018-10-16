pcrelate2 <- function(	gdsobj,
						pcs,
						scale = c('overall', 'variant', 'none'),
						ibd.probs = TRUE,
						sample.include = NULL,
						training.set = NULL,
						sample.block.size = 5000,
						snp.include = NULL,
						snp.block.size = 10000,
						maf.thresh = 0.01,
						maf.bound.method = c('truncate', 'filter'),
						small.samp.correct = FALSE,
						num.cores = 1,
						verbose = TRUE){

	# checks
	.pcrelateChecks(gdsobj = gdsobj, pcs = pcs, scale = scale, ibd.probs = ibd.probs, sample.include = sample.include, training.set = training.set, 
					maf.thresh = maf.thresh, maf.bound.method = maf.bound.method)
	
	# set up number of cores
	sys.cores <- parallel::detectCores(logical = TRUE)
	doMC::registerDoMC(cores = min(c(num.cores, sys.cores)))
	message('Using ', min(c(num.cores, sys.cores)), ' CPU cores')

	# number of sample blocks
	nsampblock <- ceiling(length(sample.include)/sample.block.size)
	# list of samples in each block
	if(nsampblock > 1){
		samp.blocks <- unname(split(sample.include, cut(1:length(sample.include), nsampblock)))
	}else{
		samp.blocks <- list(sample.include)
	}

	if(nsampblock == 1){
		message('Running analysis with ', length(sample.include), ' samples...')

		# create matrix of PCs
		V <- .createPCMatrix(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include)
		# get index of training.set
		if(!is.null(training.set)){
			tidx <- sample.include %in% training.set
			message('Betas for ', ncol(pcs), ' PC(s) will be calculated using ', sum(tidx), ' samples in training.set...')
		}else{
			tidx <- rep(TRUE, length(sample.include))
			message('Betas for ', ncol(pcs), ' PC(s) will be calculated using all ', sum(tidx), ' samples in sample.include...')
		}
		# matrix product of V
		VVtVi <- tcrossprod(V[tidx,], chol2inv(chol(crossprod(V[tidx,]))))

		### Stephanie to build in variant iterators here ###
		# number of snp blocks
		nsnpblock <- ceiling(length(snp.include)/snp.block.size)
		# list of snps in each block
		if(nsnpblock > 1){
			snp.blocks <- unname(split(snp.include, cut(1:length(snp.include), nsnpblock)))
		}else{
			snp.blocks <- list(snp.include)
		}
		message('Running PC-Relate analysis with ', length(snp.include), ' SNPs in ', nsnpblock, ' blocks...')

		# for each snp block
		matList <- foreach(k = 1:nsnpblock, .combine = .matListCombine, .inorder = FALSE, .multicombine = FALSE) %dopar% {
			# read genotype data for the block
			seqSetFilter(gdsobj, variant.id = snp.blocks[[k]])
			G <- altDosage(gdsobj)

			# calculate ISAF betas
			beta <- crossprod(G[tidx,], VVtVi)

			# calculate PC-Relate estimates
			matList.block <- .pcrelateVarBlock(	G = G, beta = beta, V = V, idx = 1:nrow(V), jdx = 1:nrow(V), scale = scale, ibd.probs = ibd.probs, 
                                    			snp.include = snp.blocks[[k]], maf.thresh = maf.thresh, maf.bound.method = maf.bound.method)
			matList.block
		}

		# take ratios to compute final estimates
		estList <- .pcrelateCalcRatio(matList = matList, scale = scale, ibd.probs = ibd.probs)

		# cast to data.tables
		dtS <- data.table(ID = rownames(estList$kin), f = 2*diag(estList$kin) - 1, nsnp = diag(estList$nsnp))
		dtB <- .estListToDT(estList, drop.lower = TRUE)


	}else if(nsampblock > 1){
		message('Splitting ', length(sample.include), ' samples in sample.include into ', nsampblock, ' blocks...')

		# calculate betas for individual specific allele frequencies
		if(!is.null(training.set)){
			message('Calculating betas for ', ncol(pcs), ' PC(s) using ', length(training.set), ' samples in training.set...')
			beta <- calcISAFBeta(gdsobj = gdsobj, pcs = pcs, sample.include = training.set, snp.include = snp.include, snp.block.size = snp.block.size)
		}else{
			message('Calculating betas for ', ncol(pcs), ' PC(s) using ', length(sample.include), ' samples in sample.include...')
			beta <- calcISAFBeta(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include, snp.include = snp.include, snp.block.size = snp.block.size)
		}
		### beta is a matrix of variants x PCs; needs to be saved, not sure of the best format or the best way to split this up ###
		
		# compute estimates for current (pair of) sample block(s)
		### this is where we would parallelize with multiple jobs ###
		dtS <- NULL
		dtB <- NULL
		for(i in 1:nsampblock){
			for(j in i:nsampblock){
				message('Computing Relatedness Estimates for sample block pair (', i, ',', j, ')')
				# compute estimates for the (pair of) sample block(s)
				out <- pcrelateSampBlock(	gdsobj = gdsobj, pcs = pcs, betaobj = beta, sample.include.block1 = samp.blocks[[i]], sample.include.block2 = samp.blocks[[j]],
											scale = scale, ibd.probs = ibd.probs, snp.include = snp.include, snp.block.size = snp.block.size, 
											maf.thresh = maf.thresh, maf.bound.method = maf.bound.method)

				# update results with this (pair of) sample block(s)
				if(i == j) dtS <- rbind(dtS, out$dtS)
				dtB <- rbind(dtB, out$dtB)
			}
		}
	}

	# order samples
	setkey(dtB, ID1, ID2)
	setkey(dtS, ID)

	### post processing after putting all samples together
	if(ibd.probs){
		# correct k2 for HW departure and use alternate k0 estimator for non-1st degree relatives
		dtB <- .fixIBDEst(dtB = dtB, dtS = dtS)
	}

	### need to rewrite the small sample correction at some point - take care of this later ###


	return(list(dtB = dtB, dtS = dtS))
}





.pcrelateChecks <- function(gdsobj, pcs, scale, ibd.probs, sample.include, training.set, maf.thresh, maf.bound.method){
	# check parameters
	if(!(scale %in% c('overall', 'variant', 'none'))) stop('scale should be one of `overall`, `variant`, or `none`')
	if(scale == 'none' & ibd.probs) stop('`ibd.probs` must be FALSE when `scale` == none')
	if(maf.thresh < 0 | maf.thresh > 0.5) stop('maf.thresh must be in [0,0.5]')
	if(!(maf.bound.method %in% c('truncate', 'filter'))) stop('maf.bound.method should be one of `truncate` or `filter`')

	# set sample filter on gdsobj
	seqSetFilter(gdsobj, sample.id = sample.include)
	# check sample.include 
	if(!all(sample.include %in% seqGetData(gdsobj, 'sample.id'))) stop('All samples in sample.include must be in the gdsobj')
	# check training.set
	if(!all(training.set %in% sample.include)) stop('All samples in training.set must be in sample.include')
	# check pcs
	if(class(pcs) != 'matrix' | is.null(rownames(pcs))) stop('pcs should be a matrix of PCs with rownames set to sample.ids')
	if(!all(sample.include %in% rownames(pcs))) stop('All samples in sample.include must be in pcs')
}




### function to match samples and create PC matrix
.createPCMatrix <- function(gdsobj, pcs, sample.include){
    # set filter to sample.include
    seqSetFilter(gdsobj, sample.id = sample.include)
    # subset and re-order pcs if needed
    V <- pcs[match(seqGetData(gdsobj, 'sample.id'), rownames(pcs)), , drop = FALSE]
    # append intercept
    V <- cbind(rep(1, nrow(V)), V)
    return(V)
}


### this function does the pcrelate estimation for a variant block
.pcrelateVarBlock <- function(G, beta, V, idx, jdx, scale, ibd.probs, snp.include, maf.thresh, maf.bound.method){

	# make sure order of G, beta, and V all line up
	if(!all.equal(rownames(G), rownames(V))) stop('G and V do not match')
	if(!all.equal(rownames(beta), colnames(G))) stop('G and beta do not match')

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
	rm(nonmiss)

	# compute kinship values
	kinList <- .pcrCalcKinOvr(G, mu, muqu, filt.idx, idx, jdx)

	if(ibd.probs){
		# compute ibd values
		ibdList <- .pcrCalcIBDOvr(G, mu, muqu, filt.idx, idx, jdx)

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

.pcrCalcIBDOvr <- function(G, mu, muqu, filt.idx, idx, jdx){
	# opposite allele frequency
	qu <- 1 - mu

	# indicator matrices of homozygotes
	Iaa <- G == 0
	IAA <- G == 2

	# dominance coded genotype matrix
	Gd <- mu
	Gd[G == 1] <- 0
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
	rm(nonmiss)

	# compute kinship values
	kinNum <- .pcrCalcKinVar(G, mu, muqu, filt.idx, idx, jdx)

	if(ibd.probs){
		# compute ibd values
		ibdList <- .pcrCalcIBDVar(G, mu, muqu, filt.idx, idx, jdx)

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

.pcrCalcIBDVar <- function(G, mu, muqu, filt.idx, idx, jdx){
	# opposite allele frequency
	qu <- 1 - mu

	# indicator matrices of homozygotes
	Iaa <- G == 0
	IAA <- G == 2

	# dominance coded genotype matrix
	Gd <- mu
	Gd[G == 1] <- 0
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
.fixIBDEst <- function(dtB, dtS){
    # correct k2 for HW departure
    dtB <- merge(dtB, dtS[,.(ID, f)], by.x = 'ID2', by.y = 'ID')
    setnames(dtB, 'f', 'f.2')
    dtB <- merge(dtB, dtS[,.(ID, f)], by.x = 'ID1', by.y = 'ID')
    setnames(dtB, 'f', 'f.1')
    dtB[, k2 := k2 - f.1*f.2][, `:=`(f.1 = NULL, f.2 = NULL)]
    
    # use alternate k0 estimator for non-1st degree relatives
    dtB[kin < 2^(-5/2), k0 := 1 - 4*kin + k2]
    
    return(dtB)
}



### exported function for computing PC betas for individual specific allele frequency calculations ###
calcISAFBeta <- function(gdsobj, pcs, sample.include, snp.include, snp.block.size){
	# checks
	# if(checks) .pcrelateChecks1(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include)

	# create matrix of PCs
	V <- .createPCMatrix(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include)
	# matrix product of V
	VVtVi <- tcrossprod(V, chol2inv(chol(crossprod(V))))
	
	### Stephanie to build in variant iterators here ###
	# number of snp blocks
	nsnpblock <- ceiling(length(snp.include)/snp.block.size)
	# list of snps in each block
	if(nsnpblock > 1){
		snp.blocks <- unname(split(snp.include, cut(1:length(snp.include), nsnpblock)))
	}else{
		snp.blocks <- list(snp.include)
	}
	message('Running analysis with ', length(snp.include), ' SNPs in ', nsnpblock, ' blocks...')

	beta <- foreach(k = 1:nsnpblock, .combine = rbind, .inorder = FALSE, .multicombine = TRUE) %dopar% {
		# read genotype data for the block
		seqSetFilter(gdsobj, variant.id = snp.blocks[[k]])
		G <- altDosage(gdsobj)

		# calculate ISAF betas
		beta.block <- crossprod(G, VVtVi)
		beta.block
	}
	### rather than returning and rbinding here, we may want to be writing the output to something

	# non-parallel version
	# beta <- NULL
	# for(k in 1:nsnpblock){
	# 	# calculate betas for PCs at each variant
	# 	beta.block <- .calcISAFBeta(gdsobj = gdsobj, VVtVi = VVtVi, snp.include = snp.blocks[[k]])
	# 	message('...SNP Block ', k, ' of ', nsnpblock, ' Completed: ', nrow(beta.block), ' SNPs')
	# 	# rbind
	# 	beta <- rbind(beta, beta.block)
	# }

	return(beta)
}


### exported function that does the pcrelate estimation for a (pair of) sample block(s)
pcrelateSampBlock <- function(gdsobj, pcs, betaobj, sample.include.block1, sample.include.block2, scale, ibd.probs, snp.include, snp.block.size, maf.thresh, maf.bound.method){

	# create (joint) PC matrix and indices
	# one sample block
	if(all(sample.include.block1 == sample.include.block2)){
		oneblock <- TRUE
		V <- .createPCMatrix(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include.block1)
		idx <- 1:nrow(V)
		jdx <- idx

	# two sample blocks
	}else{
		# check for overlapping samples
		if(any(sample.include.block1 %in% sample.include.block2)) stop('when using multiple sample blocks, each sample should be included in only one block')
		oneblock <- FALSE
		V <- .createPCMatrix(gdsobj = gdsobj, pcs = pcs, sample.include = c(sample.include.block1, sample.include.block2))
		idx <- which(rownames(V) %in% sample.include.block1)
		jdx <- which(rownames(V) %in% sample.include.block2)
	}
	### slight inefficiency above because we create the PC Matrix for samples in block1 for each block2 when we don't have to if we are running serially; 
	### but this seems more straightforward to parallelize


	# compute estimates looping over blocks of variants
	### this could be parallelized by blocks of variants in snp.include ###
	### Stephanie to build in variant iterators here ###

	# number of snp blocks
	nsnpblock <- ceiling(length(snp.include)/snp.block.size)
	# list of snps in each block
	if(nsnpblock > 1){
		snp.blocks <- unname(split(snp.include, cut(1:length(snp.include), nsnpblock)))
	}else{
		snp.blocks <- list(snp.include)
	}

	# compute estimates for each variant block; sum along the way
	matList <- foreach(k = 1:nsnpblock, .combine = .matListCombine, .inorder = FALSE, .multicombine = FALSE) %dopar% {

		# read genotype data for the block
		seqSetFilter(gdsobj, variant.id = snp.blocks[[k]])
		G <- altDosage(gdsobj)

		# load betas for the current block of variants
		beta.block <- betaobj[rownames(betaobj) %in% snp.blocks[[k]], , drop = FALSE]
		### this line of code will probably be different if we save the betas; need to load correct betas
		### rownames of beta.block needs to match colnames of G

		matList.block <- .pcrelateVarBlock(	G = G, beta = beta.block, V = V, idx = idx, jdx = jdx, scale = scale, ibd.probs = ibd.probs, 
                                    		snp.include = snp.blocks[[k]], maf.thresh = maf.thresh, maf.bound.method = maf.bound.method)
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
	# 	message('...SNP Block ', k, ' of ', nsnpblock, ' Completed: ', length(snp.blocks[[k]]), ' SNPs')
	# }

	# compute final estimates
	estList <- .pcrelateCalcRatio(matList = matList, scale = scale, ibd.probs = ibd.probs)

	# cast to data.tables
	if(oneblock){
		# self table
		dtS <- data.table(ID = rownames(estList$kin), f = 2*diag(estList$kin) - 1, nsnp = diag(estList$nsnp))
		setkey(dtS, ID)
		# between samples table
		dtB <- .estListToDT(estList, drop.lower = TRUE)
	}else{
		dtS <- NULL
		# between samples table
		dtB <- .estListToDT(estList, drop.lower = FALSE)
	}
	setkey(dtB, ID1, ID2)
	
	return(list(dtS = dtS, dtB = dtB))
}


