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

	# some parameter checks
	if(maf.thresh < 0 | maf.thresh > 0.5) stop('maf.thresh must be in [0,0.5]')
	if(scale == 'none' & ibd.probs) stop('IBD probabilities can not be calculated when `scale` == none')
	### more to add ###


	# set up number of cores
	sys.cores <- parallel::detectCores(logical = TRUE)
	doMC::registerDoMC(cores = min(c(num.cores, sys.cores)))
	message('Using ', min(c(num.cores, sys.cores)), ' CPU cores')


	# calculate betas for individual specific allele frequencies
	### this may be done externally from this function first; in this case we would load betas below ###
	if(!is.null(training.set)){
		message('Calculating betas for ', ncol(pcs), ' PC(s) using ', length(training.set), ' samples in training.set...')
		beta <- calcBetaISAF(gdsobj = gdsobj, pcs = pcs, sample.include = training.set, snp.include = snp.include, snp.block.size = snp.block.size)
	}else{
		message('Calculating betas for ', ncol(pcs), ' PC(s) using ', length(sample.include), ' samples in sample.include...')
		beta <- calcBetaISAF(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include, snp.include = snp.include, snp.block.size = snp.block.size)
	}


	# number of sample blocks
	nsampblock <- ceiling(length(sample.include)/sample.block.size)
	# list of samples in each block
	if(nsampblock > 1){
		samp.blocks <- unname(split(sample.include, cut(1:length(sample.include), nsampblock)))
		message('Splitting ', length(sample.include), ' samples in sample.include into ', nsampblock, ' blocks...')
	}else{
		samp.blocks <- list(sample.include)
	}


	# compute estimates for current (pair of) sample block(s)
	### this could be parallelized by pairs of blocks of samples ###
	### probably do this with multiple jobs? ###

	outB <- NULL
	outS <- NULL
	for(i in 1:nsampblock){
		for(j in i:nsampblock){
			message('Computing Relatedness Estimates for sample block pair (', i, ',', j, ')')
			# compute estimates for the (pair of) sample block(s)
			estList <- .pcrelateSampBlock(	gdsobj = gdsobj, beta = beta, pcs = pcs, sample.include.block1 = samp.blocks[[i]], sample.include.block2 = samp.blocks[[j]],
											scale = scale, ibd.probs = ibd.probs, snp.include = snp.include, snp.block.size = snp.block.size, 
											maf.thresh = maf.thresh, maf.bound.method = maf.bound.method)

			# update results with this (pair of) sample block(s)
			if(i == j){
				# create and rbind a data.table (self)
				outS <- rbind(outS, data.table(ID = rownames(estList$kin), f = 2*diag(estList$kin) - 1, nsnp = diag(estList$nsnp)))
				# melt and merge into a data.table (between samples)
				estDTB <- .estListToDT(estList, drop.lower = TRUE)
			}else{
				# melt and merge into a data.table (between samples)
				estDTB <- .estListToDT(estList, drop.lower = FALSE)
			}
			# rbind
			outB <- rbind(outB, estDTB)
		}
	}
	setkey(outB, ID1, ID2)
	setkey(outS, ID)

	### post processing after putting all samples together
	if(ibd.probs){
		# correct k2 for HW departure and use alternate k0 estimator for non-1st degree relatives
		outB <- .fixIBDEst(outB = outB, outS = outS)
	}

	### need to rewrite the small sample correction at some point - take care of this later ###


	return(list(outB = outB, outS = outS))
}







### function to match samples and create PC matrix
.createPCMatrix <- function(gdsobj, pcs, sample.include){
	# subset to sample.include
	seqSetFilter(gdsobj, sample.id = sample.include)
	sample.id <- seqGetData(gdsobj, 'sample.id')

	if(!all(sample.id %in% rownames(pcs))) stop('all samples must be in pcs')

	# subset and re-order pcs if needed
	V <- pcs[match(sample.id, rownames(pcs)), , drop = FALSE]
	# append intercept
	V <- cbind(rep(1, nrow(V)), V)

	return(V)
}




### these functions are for computing PC betas for individual specific allele frequency calculations ###
calcBetaISAF <- function(gdsobj, pcs, sample.include, snp.include, snp.block.size){
	# create matrix of PCs
	V <- .createPCMatrix(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include)

	# matrix product of V
	VVtVi <- tcrossprod(V, chol2inv(chol(crossprod(V))))

	### what is below could be parallelized by blocks of variants in snp.include ###
	### do we store the VVtVi matrix? it is npcs x nunrels in dimension ###

	### Stephanie to build in variant iterators here ###

	# number of snp blocks
	nsnpblock <- ceiling(length(snp.include)/snp.block.size)
	# list of snps in each block
	if(nsnpblock > 1){
		snp.blocks <- unname(split(snp.include, cut(1:length(snp.include), nsnpblock)))
	}else{
		snp.blocks <- list(snp.include)
	}

	beta <- foreach(k = 1:nsnpblock, .combine = rbind, .inorder = FALSE, .multicombine = TRUE) %dopar% {
		beta.block <- .calcBetaISAF(gdsobj = gdsobj, VVtVi = VVtVi, snp.include = snp.blocks[[k]])
		beta.block
	}

	# non-parallel version
	# beta <- NULL
	# for(k in 1:nsnpblock){
	# 	# calculate betas for PCs at each variant
	# 	beta.block <- .calcBetaISAF(gdsobj = gdsobj, VVtVi = VVtVi, snp.include = snp.blocks[[k]])
	# 	message('...SNP Block ', k, ' of ', nsnpblock, ' Completed: ', nrow(beta.block), ' SNPs')
	# 	# rbind
	# 	beta <- rbind(beta, beta.block)
	# }

	return(beta)
}

.calcBetaISAF <- function(gdsobj, VVtVi, snp.include){
	# load genotype data
	seqSetFilter(gdsobj, variant.id = snp.include)
	G <- altDosage(gdsobj)

	# compute beta estimates
	beta <- crossprod(G, VVtVi)
	### beta is a matrix of variants x PCs; probably needs to be saved, not sure of the best format ###

	return(beta)
}



### this function does the pcrelate estimation for a (pair of) sample block(s)
.pcrelateSampBlock <- function(gdsobj, beta, pcs, sample.include.block1, sample.include.block2, scale, ibd.probs, snp.include, snp.block.size, maf.thresh, maf.bound.method){

	# create (joint) PC matrix and indices

	# one sample block
	if(all(sample.include.block1 == sample.include.block2)){
		V <- .createPCMatrix(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include.block1)
		idx <- 1:nrow(V)
		jdx <- idx

	# two sample blocks
	}else{
		# check for overlapping samples
		if(any(sample.include.block1 %in% sample.include.block2)) stop('when using multiple sample blocks, each sample should be included in only one block')
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
		matList.block <- .pcrelateVarBlock(	gdsobj = gdsobj, beta = beta, V = V, idx = idx, jdx = jdx, scale = scale, ibd.probs = ibd.probs, 
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

.matListCombine <- function(...){
    mapply(FUN = "+", ..., SIMPLIFY = FALSE)
}



### this function does the pcrelate estimation for a variant block
.pcrelateVarBlock <- function(gdsobj, beta, V, idx, jdx, scale, ibd.probs, snp.include, maf.thresh, maf.bound.method){

	# load genotype data
	seqSetFilter(gdsobj, variant.id = snp.include, sample.id = rownames(V))
	G <- altDosage(gdsobj)
	### rownames(V) needs to match rownames(G) - it should based on the .createPCMatrix() function

	# load betas for the current block of variants (colnames(G))
	beta.block <- beta[match(colnames(G), rownames(beta)), , drop = FALSE]
	### this line of code will probably be different if we save the betas; need to load correct betas
	### rownames of beta.block needs to match colnames of G

	# estimate individual specific allele frequencies
	mu <- .estISAF(V = V, beta = beta.block, bound.thresh = maf.thresh, bound.method = maf.bound.method)

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



### these functions are for estimating individual specific allele frequencies
.estISAF <- function(V, beta, bound.thresh, bound.method){
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


### functions for post processing
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

.fixIBDEst <- function(outB, outS){
    # correct k2 for HW departure
    outB <- merge(outB, outS[,.(ID, f)], by.x = 'ID2', by.y = 'ID')
    setnames(outB, 'f', 'f.2')
    outB <- merge(outB, outS[,.(ID, f)], by.x = 'ID1', by.y = 'ID')
    setnames(outB, 'f', 'f.1')
    outB[, k2 := k2 - f.1*f.2][, `:=`(f.1 = NULL, f.2 = NULL)]
    
    # use alternate k0 estimator for non-1st degree relatives
    outB[kin < 2^(-5/2), k0 := 1 - 4*kin + k2]
    
    return(outB)
}