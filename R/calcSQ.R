

.computeSigmaQuantities <- function(varComp, covMatList, group.idx = NULL, vmu = NULL, gmuinv = NULL){
    m <- length(covMatList)
    n <- nrow(covMatList[[1]])

###    Sigma <- Vre <- Reduce("+", mapply("*", covMatList, varComp[1:m], SIMPLIFY=FALSE))
    Vre <- Reduce("+", mapply("*", covMatList, varComp[1:m], SIMPLIFY=FALSE))
    
    if (is.null(vmu)){ ## this means the family is "gaussian"
        if (is.null(group.idx)){
            Sigma <- Diagonal(x=rep(varComp[m+1],n)) + Vre
        } else{
            g <- length(group.idx)
        
###        diagV <- rep(0,nrow(covMatList[[1]]))
###        for(i in 1:g){
###            diagV[group.idx[[i]]] <- varComp[m+i]
###        }
###        diag(Sigma) <- diag(Sigma) + diagV

            mylevels <- rep(NA, n)
            for(i in 1:g){
                mylevels[group.idx[[i]]] <- i # setting up vector of indicators; who is in which group
            }
            Sigma <- Diagonal(x=(varComp[m+1:g])[mylevels] ) + Vre
        }

    ### if non-gaussian family:
    } else {
###        Sigma <- Sigma + diag(as.vector(vmu)/as.vector(gmuinv)^2)
        Sigma <- Diagonal(x=as.vector(vmu)/as.vector(gmuinv)^2) + Vre
    }
    
    # cholesky decomposition
    cholSigma <- chol(Sigma)
    # inverse
    Sigma.inv <- chol2inv(cholSigma)

    return(list(cholSigma = cholSigma, Sigma.inv = Sigma.inv, Vre = Vre))

}

