
calcBetaISAF <- function(gdsobj, pcs, training.set, snp.include){

	V <- .createPCMatrix(gdsobj = gdsobj, pcs = pcs, sample.include = training.set)

	# matrix product of V
	VVtVi <- tcrossprod(V, chol2inv(chol(crossprod(V))))
	### do we store this matrix? it is npcs x nunrels in dimension ###

	### this could be parallelized by blocks of variants in snp.include ###
	.calcBetaISAF(gdsobj, snp.include, VVtVi)

}

.calcBetaISAF <- function(gdsobj, snp.include, VVtVi){

	# load genotype data
	seqSetFilter(gdsobj, variant.sel = snp.include)
	G <- altDosage(gdsobj)

	# compute beta estimates
	beta <- crossprod(G, VVtVi)
	### beta is a matrix of variants x PCs; needs to be saved, not sure of the best format ###
}

.createPCMatrix <- function(gdsobj, pcs, sample.include){
	# subset to sample.include
	seqSetFilter(gdsobj, sample.id = sample.include)
	sample.id <- seqGetData(gdsobj, 'sample.id')

	if(!all(sample.id %in% rownames(pcs))){
		stop('all samples must be in pcs')
	}

	# subset and re-order pcs if needed; append intercept
	V <- pcs[match(sample.id, rownames(pcs)),]
	V <- cbind(rep(1, nrow(V)), V)

	return(V)
}



### kinship for one block of samples and one block of variants ###


V <- .createPCMatrix(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include)


### this could be parallelized by blocks of variants in snp.include ###


# load genotype data
seqSetFilter(gdsobj, variant.sel = snp.include, sample.id = rownames(V))
G <- altDosage(gdsobj)

# index of missing values
filt.idx <- which(is.na(G))

# get ISAF estimates (i.e. 0.5*fitted values)
mu <- 0.5*tcrossprod(V, beta)
# fix boundary cases
mu <- apply(mu, 2, .muBoundaryFn, thresh = 0.025)

if(scale != 'none'){
	# calculate mu(1-mu)
	muqu <- mu*(1 - mu)
}


if(scale = 'overall'){
	kinList <- .calcKinIndOvr(G, mu, muqu, filt.idx)
}else if(scale == 'variant'){
	kin <- .calcKinIndVar(G, mu, muqu, filt.idx)
}else if(scale == 'none'){
	kin <- .calcKinNone(G, mu, filt.idx)
}

if(ibd.probs){
	if(scale == 'overall'){
		ibdList <- .calcIBDIndOvr(G, mu, muqu, filt.idx)
	}else if(scale == 'variant'){
		ibdList <- .calcIBDIndVar(G, mu, muqu, filt.idx)
	}
}



.calcKinIndOvr <- function(G, mu, muqu, filt.idx){
	# residuals
	R <- G - 2*mu
	# set missing values to 0 (i.e. no contribution from that variant)
	R[filt.idx] <- 0
	# numerator (crossprod of residuals)
	kinNum <- tcrossprod(R)

	# set missing value to 0 (i.e. no contribution from that variant)
	muqu[filt.idx] <- 0
	# denominator (needs to be multipled by 2 when all put together)
	kinDen <- tcrossprod(sqrt(muqu))

	return(list(kinNum = kinNum, kinDen = kinDen))
}

.calcKinIndVar <- function(G, mu, muqu, filt.idx){
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

.calcIBDIndOvr <- function(G, mu, muqu, filt.idx){
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

.calcIBDIndVar <- function(G, mu, muqu, filt.idx){
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


# function to fix boundary issues with ISAFs
.muBoundaryFn <- function(x, thresh){
    x[x < thresh] <- thresh
    x[x > (1 - thresh)] <- (1 - thresh)
    return(x)
}









##### I'd like to just get rid of the population frequency based calculations #####



### freq.type == 'population' ###
p <- 0.5*colMeans(G[rownames(G) %in% training.set,], na.rm = TRUE)

if(scale != 'none'){
	# calculate p(1-p)
	pq <- p*(1-p)
}


.calcKinPopOvr <- function(G, p, pq){
	# residuals
	R <- sweep(G, MARGIN = 2, STATS = 2*p, FUN = '-')
	# set missing values to 0 (i.e. no contribution from that variant)
	R[filt.idx] <- 0
	# numerator (crossprod of residuals)
	kinNum <- tcrossprod(R)

	# create matrix for denominator computation
	pqMat <- matrix(sqrt(pq), nrow = nrow(G), ncol = ncol(G), byrow = TRUE)
	# set missing value to 0 (i.e. no contribution from that variant)
	pqMat[filt.idx] <- 0
	# denominator (needs to be multipled by 2 when all put together)
	kinDen <- tcrossprod(pqMat)

	return(list(kinNum = kinNum, kinDen = kinDen))
}

.calcKinPopVar <- function(G, p, pq){
	# residuals
	R <- sweep(G, MARGIN = 2, STATS = 2*p, FUN = '-')
	# scaled residuals
	R <- sweep(R, MARGIN = 2, STATS = sqrt(pq), FUN = '/')
	# set missing values to 0 (i.e. no contribution from that variant)
	R[filt.idx] <- 0
	# crossprod of scaled residuals
	kin <- tcrossprod(R)

	return(kin)
}

.calcIBDPopOvr <- function(G, p, pq){
	# indicator matrices of homozygotes
	Iaa <- G == 0
	IAA <- G == 2

	# dominance coded genotype matrix
	Gd <- matrix(p, nrow = nrow(G), ncol = ncol(G), byrow = TRUE)
	Gd[G == 1] <- 0
	Gd[IAA] <- 1 - Gd[IAA]
	Gd <- sweep(Gd, MARGIN = 2, STATS = pq, FUN = '-')

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







