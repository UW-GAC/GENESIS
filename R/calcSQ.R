# .computeSigmaQuantities <- function(varComp, covMatList, group.idx = NULL, vmu = NULL, gmuinv = NULL){
.computeSigmaQuantities <- function(varComp, covMatList, group.idx = NULL, diagV = NULL){
    m <- length(covMatList)
    n <- nrow(covMatList[[1]])

    # contribution to Sigma from random effects
    Vre <- Reduce("+", mapply("*", covMatList, varComp[1:m], SIMPLIFY=FALSE))

    if (is.null(diagV)){
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
        }
    }

    # NEED TO MULTIPLY diagV BY varComp[m+1] WHEN GAUSSIAN

    # } else {
    #     # non-gaussian family
    #     diagV <- as.vector(vmu)/as.vector(gmuinv)^2
    # }

    # construct Sigma
    Sigma <- Vre
    diag(Sigma) <- diag(Sigma) + diagV
    # cholesky decomposition
    cholSigma <- chol(Sigma)
    # inverse
    Sigma.inv <- chol2inv(cholSigma)

    return(list(Sigma.inv = Sigma.inv, Vre = Vre, W = 1/diagV, cholSigma.diag = diag(cholSigma)))
    # W is the diagonal matrix needed for fast.score.SE
}
