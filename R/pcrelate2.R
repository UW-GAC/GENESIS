
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
	seqSetFilter(gdsobj, variant.sel = snp.include)
	G <- altDosage(gdsobj)

	# get ISAF estimates (i.e. 0.5*fitted values)
	mu <- 0.5*tcrossprod(V, beta)
	# fix boundary cases
	mu <- apply(mu, 2, .muBoundaryFn, thresh = 0.025)

	# numerator of kinship matrix (crossprod of residuals)
	kinNum <- tcrossprod(G - 2*mu)

	# calculate mu(1-mu)
	muqu <- mu*(1-mu)
	# denominator of kinship marix (needs to be multipled by 2 when all put together)
	kinDen <- tcrossprod(sqrt(muqu))

	if(ibd.probs){
		# k0
		# indicator matrices of homozygotes
		Iaa <- G == 0
		IAA <- G == 2
    	k0Num <- tcrossprod(IAA, Iaa) + tcrossprod(Iaa, IAA)    	
    	rm(Iaa)
    	k0Den <- tcrossprod(mu^2, (1-mu)^2) + tcrossprod((1-mu)^2, mu^2)

		# k2
		# dominance coded genotype matrix
		Gd <- mu
    	Gd[G == 1] <- 0
    	Gd[IAA] <- 1 - mu[IAA] 
    	Gd <- Gd - muqu
    	k2Num <- tcrossprod(Gd)
    	k2Den <- tcrossprod(muqu)

    	
	}




# function to fix boundary issues with ISAFs
.muBoundaryFn <- function(x, thresh){
    x[x < thresh] <- thresh
    x[x > (1 - thresh)] <- (1 - thresh)
    return(x)
}
