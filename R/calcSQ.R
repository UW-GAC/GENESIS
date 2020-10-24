.computeSigmaQuantities <- function(varComp, covMatList, group.idx = NULL, vmu = NULL, gmuinv = NULL, nsamp, ntraits=1){
    m <- length(covMatList)
#    n <- nrow(covMatList[[1]])
    n <- nsamp

    ## contribution to Sigma from random effects
    if(ntraits==1){
    Vre <- Reduce("+", mapply("*", covMatList, varComp[1:m], SIMPLIFY=FALSE))
    }else{
        kinMatSmall <- covMatList[[1]][1:nsamp,1:nsamp]
        Vre <- matrix(NA,nrow=ntraits*nsamp,ncol=ntraits*nsamp)
        Vre[1:nsamp,1:nsamp] <- matrix(varComp[3]*kinMatSmall)
        Vre[1:nsamp,(nsamp+1):(ntraits*nsamp)] <- matrix(varComp[5]*kinMatSmall)
        Vre[(nsamp+1):(ntraits*nsamp),1:nsamp] <- matrix(varComp[5]*kinMatSmall)
        Vre[(nsamp+1):((ntraits*nsamp)),(nsamp+1):((ntraits*nsamp))] <- matrix(varComp[4]*kinMatSmall)
    }
    
    if (is.null(vmu)){
        # gaussian family
        # contribution to Sigma from residual variance
        if (is.null(group.idx)){
            diagV <- rep(varComp[m+1],n)
        } else{
            g <- length(group.idx)
            diagV <- rep(NA, n)
            for(i in 1:g){
                diagV[group.idx[[i]]] <- varComp[m+i]
            }
            
            diagV[1:nsamp] <- varComp[1]
            diagV[(nsamp+1):(nsamp*ntraits)] <- varComp[2]
            
            # mylevels <- rep(NA, n)
            # for(i in 1:g){
            #     mylevels[group.idx[[i]]] <- i # setting up vector of indicators; who is in which group
            # }
            # diagV <- (varComp[(m+1):(m+g)])[mylevels]
        }

    } else {
        # non-gaussian family
        diagV <- as.vector(vmu)/as.vector(gmuinv)^2
    }

    # construct Sigma
    Sigma <- Vre
    diag(Sigma) <- diag(Sigma) + diagV   
    
    # check for PD, if not, estimate to nearest PD matrix
    # we need Sigma to be PD to take chol decomposition
    if(!is.positive.definite(Sigma)){
        Sigma <- make.positive.definite(Sigma)
    }
    # cholesky decomposition
    cholSigma <- chol(Sigma)
    # inverse
    Sigma.inv <- chol2inv(cholSigma)

    return(list(Sigma.inv = Sigma.inv, Vre = Vre, diagV = diagV, cholSigma.diag = diag(cholSigma)))

}
