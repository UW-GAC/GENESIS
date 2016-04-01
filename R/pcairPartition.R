pcairPartition <- function(kinMat, kin.thresh = 2^(-11/2), divMat=NULL, div.thresh = -2^(-11/2), unrel.set=NULL){	

    # check that kinMat has colummn names
    if(is.null(colnames(kinMat))){
    	stop("colnames and rownames of kinMat must be individual IDs")
    }
	# create vector of IDs
	IDs <- colnames(kinMat)
	
	# check that row and column names match
	if(!all(IDs == rownames(kinMat))){
		stop("colnames and rownames of kinMat do not match")
	}

	# set diagonal to 0
	diag(kinMat) <- 0

	# check that unrel.set is included in kinMat
	if(!is.null(unrel.set)){
		if(!all(unrel.set %in% IDs)){
			stop("All of the samples in unrel.set must be in kinMat")
		}
	}
	
	# if divMat unspecified
	if(is.null(divMat)){
		message("divMat not specified; using kinMat for kinship and divergence measures...")
		divMat <- kinMat
	# if divMat specified
	}else{
		# check that colnames and rownames match kinMat
		if(!all(colnames(divMat) == rownames(divMat))){
			stop("colnames and rownames of divMat do not match")
		}
		if(!all(IDs == colnames(divMat))){
			stop("colnames of kinMat and divMat do not match")
		}
	}


	# indicator for divergent pairs
	indDIV <- (divMat < -abs(div.thresh)) & (kinMat < kin.thresh)
	rm(divMat)
	# number of divergent pairs
	ndiv <- colSums(indDIV)
	rm(indDIV)

	# "total kinship" values
	kinMat[kinMat < kin.thresh] <- 0
	kinsum <- colSums(kinMat)
	# indicator for relatives
	indKIN <- kinMat != 0
	rm(kinMat)
	
	
	ridx <- NULL  # index of samples in related set
	uidx <- NULL  # index of samples in unrelated set
	if(!is.null(unrel.set)){
		# index of user specified unrelated set
		unrel.idx <- which(IDs %in% unrel.set)
		# index of relatives of the unrelated set
		rel.idx <- as.vector(which(colSums(indKIN[unrel.idx,]) > 0))
		# remove individuals in unrel.idx from rel.idx
		rel.idx <- rel.idx[!(rel.idx %in% unrel.idx)]

		# append
		ridx <- append(ridx, rel.idx)
		uidx <- append(uidx, unrel.idx)
		# subset indKIN
		indKIN <- indKIN[-c(rel.idx, unrel.idx),-c(rel.idx, unrel.idx)]

		rm(unrel.idx); rm(rel.idx)
	}

    # number of relatives	
	numrel <- colSums(indKIN)

	# identify individuals with no relatives
	unrel.idx <- which(IDs %in% colnames(indKIN)[which(numrel == 0)])	
	uidx <- append(uidx, unrel.idx)
	
	# subset indKIN to individuals with relatives
	keep.idx <- which(numrel > 0)
	indKIN <- indKIN[keep.idx, keep.idx]

	# find connected components
	grph <- as(indKIN, "graphNEL")
	conn <- connComp(grph)

	# IDs of individuals left in indKIN
	subIDs <- colnames(indKIN)
	

	# loop through components
	for(i in 1:length(conn)){		

		# indices for IDs
		ID.idx <- which(IDs %in% conn[[i]])
		subID.idx <- which(subIDs %in% conn[[i]])

		# subset indKIN
		indKINsub <- indKIN[subID.idx, subID.idx]
		
		# number of relatives
		numrel <- colSums(indKINsub)
		rel.idx <- NULL

		# partition the component		
		while(max(numrel) > 0){
			idx <- which(numrel == max(numrel))
			if(length(idx) > 1){
				ndiv2 <- ndiv[ID.idx][idx]
				didx <- which(ndiv2 == min(ndiv2))
				if(length(didx) == 1){
					idx <- idx[didx]
				}else{
					kinsum2 <- kinsum[ID.idx][idx[didx]]
					kidx <- which(kinsum2 == min(kinsum2))[1]
					idx <- idx[didx[kidx]]
				}
			}
			# update
			rel.idx <- append(rel.idx, ID.idx[idx])
			indKINsub[idx,] <- FALSE; indKINsub[,idx] <- FALSE
			numrel <- colSums(indKINsub)
		}

		# unrelateds from this component
		unrel.idx <- ID.idx[!(ID.idx %in% rel.idx)]

		# append
		ridx <- append(ridx, rel.idx)
		uidx <- append(uidx, unrel.idx)
	}


	# collect results
	if(!is.null(ridx)){
		ridx <- as.numeric(sort(ridx))
		rels <- IDs[ridx]
	}else{
		rels <- NULL
	}
	uidx <- as.numeric(sort(uidx))
	unrels <- IDs[uidx]
		
	# return results
	return(list(rels = rels, unrels=unrels))
}
