
calcBetaISAF <- function(gdsobj, pcs, sample.include, snp.include){
	# create matrix of PCs
	V <- .createPCMatrix(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include)

	# matrix product of V
	VVtVi <- tcrossprod(V, chol2inv(chol(crossprod(V))))
	### do we store this matrix? it is npcs x nunrels in dimension ###

	### this could be parallelized by blocks of variants in snp.include ###
	# calculate betas for PCs at each variant
	.calcBetaISAF(gdsobj, snp.include, VVtVi)
}

.createPCMatrix <- function(gdsobj, pcs, sample.include){
	# subset to sample.include
	seqSetFilter(gdsobj, sample.id = sample.include)
	sample.id <- seqGetData(gdsobj, 'sample.id')

	if(!all(sample.id %in% rownames(pcs))){
		stop('all samples must be in pcs')
	}

	# subset and re-order pcs if needed
	V <- pcs[match(sample.id, rownames(pcs)), , drop = FALSE]
	# append intercept
	V <- cbind(rep(1, nrow(V)), V)

	return(V)
}

.calcBetaISAF <- function(gdsobj, VVtVi, snp.include){
	# load genotype data
	seqSetFilter(gdsobj, variant.sel = snp.include)
	G <- altDosage(gdsobj)

	# compute beta estimates
	beta <- crossprod(G, VVtVi)
	### beta is a matrix of variants x PCs; needs to be saved, not sure of the best format ###
}





### kinship for one block of samples and one block of variants ###


V <- .createPCMatrix(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include)


### this could be parallelized by blocks of variants in snp.include ###


# load genotype data
seqSetFilter(gdsobj, variant.sel = snp.include, sample.id = rownames(V))
G <- altDosage(gdsobj)
### rownames(V) needs to match rownames(G) - it should based on the .createPCMatrix() function

# load betas for the current block of variants (colnames(G))
beta.block <- # something
### rownames of beta.block needs to match colnames of G

# estimate individual specific allele frequencies
mu <- .estISAF(V, beta.block, bound.method = bound.method, bound.thresh = bound.thresh)

.estISAF <- function(V, beta, bound.method = c('truncate','filter'), bound.thresh){
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


if(scale != 'none'){
	# calculate mu(1-mu)
	muqu <- mu*(1 - mu)
}


if(scale == 'variant'){
	# compute number of overlapping snps by pair
	nonmiss <- !(is.na(G) | is.na(mu))
	nsnp <- crossprod(nonmiss)
	# index of missing values
	filt.idx <- which(!nonmiss)
	rm(nonmiss)
}else{
	# index of missing values
	filt.idx <- which(is.na(G) | is.na(mu))
}


filt.idx <- which(is.na(G) | is.na(mu))





# calculate kinship values for the block of variants
if(scale = 'overall'){
	kinList <- .calcKinOvr(G, mu, muqu, filt.idx)
}else if(scale == 'variant'){
	kin <- .calcKinVar(G, mu, muqu, filt.idx)
}else if(scale == 'none'){
	kin <- .calcKinNone(G, mu, filt.idx)
}

# calculate ibd values for the block of variants
if(ibd.probs){
	if(scale == 'overall'){
		ibdList <- .calcIBDOvr(G, mu, muqu, filt.idx)
	}else if(scale == 'variant'){
		ibdList <- .calcIBDVar(G, mu, muqu, filt.idx)
	}
}


### this needs to be done outside of the parallelization; all the matrices need to be summed first ###

# compute final estimates
if(scale == 'overall'){
	kin <- kinList$kinNum/(4*kinList$kinDen)
}else if(scale == 'variant'){
	kin <- kin/(4*nsnp)
}

if(ibd.probs){
	if(scale == 'overall'){
		k2 <- ibdList$k2Num/ibdList$k2Den
		k0 <- ibdList$k0Num/ibdList$k0Den
	}else if(scale == 'variant'){
		k2 <- ibdList$k2/nsnp
		k0 <- ibdList$k0/nsnp
	}

	# correct k2 for HW departure
	f <- 2*diag(kin) - 1
	fprod <- tcrossprod(f)
	k2 <- k2 - fprod

	# use alternate k0 estimator for non-1st degree relatives
	k0.idx <- which(kin < 2^(-5/2))
	k0[k0.idx] <- 1 - 4*kin[k0.idx] + k2[k0.idx]
}



.calcKinOvr <- function(G, mu, muqu, filt.idx){
	# residuals
	R <- G - 2*mu	
	# set missing values to 0 (i.e. no contribution from that variant)
	R[filt.idx] <- 0
	muqu[filt.idx] <- 0
	# numerator (crossprod of residuals)
	kinNum <- tcrossprod(R)
	# denominator (needs to be multipled by 2 when all put together)
	kinDen <- tcrossprod(sqrt(muqu))

	return(list(kinNum = kinNum, kinDen = kinDen))
}

.calcKinVar <- function(G, mu, muqu, filt.idx){
	# scaled residuals
	R <- (G - 2*mu)/sqrt(muqu)
	# set missing values to 0 (i.e. no contribution from that variant)
	R[filt.idx] <- 0
	# crossprod of scaled residuals
	kin <- tcrossprod(R)

	return(kin)
}

.calcKinNone <- function(G, mu, filt.idx){
	# residuals
	R <- G - 2*mu
	# set missing values to 0 (i.e. no contribution from that variant)
	R[filt.idx] <- 0
	# crossprod of scaled residuals
	kin <- tcrossprod(R)

	return(kin)
}

.calcIBDOvr <- function(G, mu, muqu, filt.idx){
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

.calcIBDVar <- function(G, mu, muqu, filt.idx){
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
	k0 <- tcrossprod(IAA, Iaa) + tcrossprod(Iaa, IAA)

	# k2
	k2 <- tcrossprod(Gd)

	return(list(k0 = k0, k2 = k2))
}





