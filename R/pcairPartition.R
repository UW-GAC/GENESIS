pcairPartition <- function(kinMat, kin.thresh = 0.025, divMat=NULL, div.thresh = -0.025, unrel.set=NULL){

    # check that kinMat has colummn names
    if(!is.null(kinMat)){
        if(is.null(colnames(kinMat))){
            stop("colnames and rownames of kinMat must be individual IDs")
        }
    }
	# create vector of IDs
	IDs <- colnames(kinMat)
	
	# check that row and column names match
	if(!all(IDs == rownames(kinMat))){
		stop("colnames and rownames of kinMat do not match")
	}
		
	# set diagonal to 0
	diag(kinMat) <- 0
	
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
	
	# check that unrel.set is included in kinMat
	if(!is.null(unrel.set)){
		if(!all(unrel.set %in% IDs)){
			stop("All of the samples in unrel.set must be in kinMat")
		}
	}
	
	
	# "total kinship" values
	kinMat[kinMat < kin.thresh] <- 0
	kinsum <- rowSums(kinMat)
	# indicator for relatives
	indKIN <- kinMat != 0
	
	# number of divergent pairs
	# zero out divergence measures for relatives
	divMat[indKIN] <- 0
	ndiv <- apply(divMat,1,function(x){ sum(x < -abs(div.thresh)) })
	
	
	ridx <- NULL  # index of samples in related set
	if(!is.null(unrel.set)){
		unrel.idx <- which(IDs %in% unrel.set)
		for(i in 1:length(unrel.idx)){
			# identify relatives of individual set to the unrelated set
			rel.idx <- which(indKIN[unrel.idx[i],])
			# don't count relatives that are also set to the unrelated set
			rel.idx <- rel.idx[!is.element(rel.idx,unrel.idx)]
			# add these relatives to the realated set
			ridx <- append(ridx, rel.idx)
			indKIN[rel.idx,] <- FALSE; indKIN[,rel.idx] <- FALSE
		}
		# ignore relationships between individuals in unrel.set
		indKIN[unrel.idx,unrel.idx] <- FALSE
	}

    # number of relatives	
	numrel <- rowSums(indKIN)
	while(max(numrel) > 0){
		idx <- which(numrel == max(numrel))
		if(length(idx) > 1){
			ndiv2 <- ndiv[idx]
			didx <- which(ndiv2 == min(ndiv2))
			if(length(didx) == 1){
				idx <- idx[didx]
			}else{
				kinsum2 <- kinsum[idx[didx]]
				kidx <- which(kinsum2 == min(kinsum2))[1]
				idx <- idx[didx[kidx]]
			}
		}
		ridx <- append(ridx, idx)
		indKIN[idx,] <- FALSE; indKIN[,idx] <- FALSE;
		numrel <- rowSums(indKIN)
	}
	
	
	# collect results
	if(!is.null(ridx)){
		ridx <- as.numeric(sort(ridx))
		rels <- IDs[ridx]
		unrels <- IDs[!(IDs %in% rels)]
	}else{
		rels <- NULL
		unrels <- IDs
	}	
	
	# return results
	return(list(rels = rels, unrels=unrels))
}
