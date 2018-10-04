pcrelate2 <- function(	gdsobj,
						pcs,
						scale = c('overall', 'variant', 'none'),
						ibd.probs = TRUE,
						sample.include = NULL,
						training.set = NULL,
						snp.include = NULL,
						maf.thresh = 0.01,
						maf.bound.method = c('truncate', 'filter'),
						small.samp.correct = FALSE,
						verbose = TRUE){

	# some parameter checks
	if(maf.thresh < 0 | maf.thresh > 0.5) stop('maf.thresh must be in [0,0.5]')
	if(scale == 'none' & ibd.probs) stop('IBD probabilities can not be calculated when `scale` == none')
	### more to add ###


	# calculate betas for individual specific allele frequencies
	### this may be done externally from this function first; in this case we would load betas below ###
	if(!is.null(training.set)){
		beta <- calcBetaISAF(gdsobj = gdsobj, pcs = pcs, sample.include = training.set, snp.include = snp.include)
	}else{
		beta <- calcBetaISAF(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include, snp.include = snp.include)
	}


	# get PC matrix for current block of samples
	V <- .createPCMatrix(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include)
	

	# compute estimates for current blocks of samples

	# compute estimates looping over blocks of variants
	### this could be parallelized by blocks of variants in snp.include ###
	### Stephanie to build in variant iterators here ###
	for(i in 1:nsnpblock){
		matList <- .pcrelateVarBlock(	gdsobj = gdsobj, beta = beta, V = V, scale = scale, ibd.probs = ibd.probs, 
										snp.include = snp.block, maf.thresh = maf.thresh, maf.bound.method = maf.bound.method)
	}


	### this needs to be done outside of the parallelization; all the matrices need to be summed first ###

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

		# correct k2 for HW departure
		f <- 2*diag(kin) - 1
		fprod <- tcrossprod(f)
		k2 <- k2 - fprod

		# use alternate k0 estimator for non-1st degree relatives
		k0.idx <- which(kin < 2^(-5/2))
		k0[k0.idx] <- 1 - 4*kin[k0.idx] + k2[k0.idx]
	}

	### need to rewrite the small sample correction at some point - take care of this later ###


	return(list(kin = kin, k0 = k0, k2 = k2))
}





### these functions are for computing PC betas for individual specific allele frequency calculations ###


calcBetaISAF <- function(gdsobj, pcs, sample.include, snp.include){
	# create matrix of PCs
	V <- .createPCMatrix(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include)

	# matrix product of V
	VVtVi <- tcrossprod(V, chol2inv(chol(crossprod(V))))
	### do we store this matrix? it is npcs x nunrels in dimension ###

	### this could be parallelized by blocks of variants in snp.include ###
	# calculate betas for PCs at each variant
	beta <- .calcBetaISAF(gdsobj = gdsobj, VVtVi = VVtVi, snp.include = snp.include)

	return(beta)
}


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

.calcBetaISAF <- function(gdsobj, VVtVi, snp.include){
	# load genotype data
	seqSetFilter(gdsobj, variant.id = snp.include)
	G <- altDosage(gdsobj)

	# compute beta estimates
	beta <- crossprod(G, VVtVi)
	### beta is a matrix of variants x PCs; needs to be saved, not sure of the best format ###

	return(beta)
}


### this function does the pcrelate estimation for a variant block

.pcrelateVarBlock <- function(gdsobj, beta, V, scale, ibd.probs, snp.include, maf.thresh, maf.bound.method){

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
		matList <- .pcrCalcOvr(G = G, mu = mu, ibd.probs = ibd.probs)
	}else if(scale == 'variant'){
		matList <- .pcrCalcVar(G = G, mu = mu, ibd.probs = ibd.probs)
	}else if(scale == 'none'){
		matList <- .pcrCalcNone(G = G, mu = mu)
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

.pcrCalcOvr <- function(G, mu, ibd.probs){
	# compute mu(1 - mu)
	muqu <- mu*(1 - mu)

	# index of missing values
	filt.idx <- which(is.na(G) | is.na(mu))

	# compute kinship values
	kinList <- .pcrCalcKinOvr(G, mu, muqu, filt.idx)

	if(ibd.probs){
		# compute ibd values
		ibdList <- .pcrCalcIBDOvr(G, mu, muqu, filt.idx)

		return(list(kinNum = kinList$kinNum, 
					kinDen = kinList$kinDen,
					k0Num = ibdList$k0Num,
					k0Den = ibdList$k0Den,
					k2Num = ibdList$k2Num,
					k2Den = ibdList$k2Den))

	}else{
		return(kinList)
	}
}

.pcrCalcKinOvr <- function(G, mu, muqu, filt.idx){
	# residuals
	R <- G - 2*mu

	# set missing values to 0 (i.e. no contribution from that variant)
	R[filt.idx] <- 0
	muqu[filt.idx] <- 0

	# numerator (crossprod of residuals)
	kinNum <- tcrossprod(R)
	# denominator
	kinDen <- tcrossprod(sqrt(muqu))

	return(list(kinNum = kinNum, kinDen = kinDen))
}

.pcrCalcIBDOvr <- function(G, mu, muqu, filt.idx){
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
	k0Num <- tcrossprod(IAA, Iaa) + tcrossprod(Iaa, IAA)
	k0Den <- tcrossprod(mu^2, qu^2) + tcrossprod(qu^2, mu^2)

	# k2
	k2Num <- tcrossprod(Gd)
	k2Den <- tcrossprod(muqu)

	return(list(k0Num = k0Num, k0Den = k0Den, k2Num = k2Num, k2Den = k2Den))
}


.pcrCalcVar <- function(G, mu, ibd.probs){
	# compute mu(1 - mu)
	muqu <- mu*(1 - mu)

	# compute number of observed snps by pair
	nonmiss <- !(is.na(G) | is.na(mu))
	nsnp <- tcrossprod(nonmiss)

	# index of missing values
	filt.idx <- which(!nonmiss)
	rm(nonmiss)

	# compute kinship values
	kinNum <- .pcrCalcKinVar(G, mu, muqu, filt.idx)

	if(ibd.probs){
		# compute ibd values
		ibdList <- .pcrCalcIBDVar(G, mu, muqu, filt.idx)

		return(list(kinNum = kinNum,
					k0Num = ibdList$k0Num,
					k2Num = ibdList$k2Num,
					nsnp = nsnp))

	}else{
		return(list(kinNum = kinNum, nsnp = nsnp))
	}
}

.pcrCalcKinVar <- function(G, mu, muqu, filt.idx){
	# scaled residuals
	R <- (G - 2*mu)/sqrt(muqu)

	# set missing values to 0 (i.e. no contribution from that variant)
	R[filt.idx] <- 0

	# numerator (crossprod of scaled residuals)
	kinNum <- tcrossprod(R)
	
	return(kinNum)
}

.pcrCalcIBDVar <- function(G, mu, muqu, filt.idx){
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
	k0Num <- tcrossprod(IAA, Iaa) + tcrossprod(Iaa, IAA)

	# k2
	k2Num <- tcrossprod(Gd)

	return(list(k0Num = k0Num, k2Num = k2Num))
}


.pcrCalcNone <- function(G, mu){
	# index of missing values
	filt.idx <- which(is.na(G) | is.na(mu))

	# compute kinship values
	kin <- .pcrCalcKinNone(G, mu, filt.idx)
	return(kin)
}

.pcrCalcKinNone <- function(G, mu, filt.idx){
	# residuals
	R <- G - 2*mu

	# set missing values to 0 (i.e. no contribution from that variant)
	filt.idx <- which(is.na(G) | is.na(mu))	
	R[filt.idx] <- 0

	# crossprod of scaled residuals
	kin <- tcrossprod(R)

	return(kin)
}



