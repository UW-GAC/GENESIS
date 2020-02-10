.calcLikelihoodQuantities <- function(Y, X, Sigma.inv, cholSigma.diag){

    if (is(Sigma.inv, "Matrix")) X <- Matrix(X)
    n <- length(Y)
    k <- ncol(X)
    
    ### Calulate the generalized least squares estimate
    Sigma.inv_X <- crossprod(Sigma.inv, X)
    Xt_Sigma.inv_X <- crossprod(X, Sigma.inv_X)
    # fix issue with not recognizing the matrix as symmetric
    Xt_Sigma.inv_X <- (Xt_Sigma.inv_X + t(Xt_Sigma.inv_X))/2
    chol.Xt_Sigma.inv_X <- chol(Xt_Sigma.inv_X)
    Xt_Sigma.inv_X.inv <- chol2inv(chol.Xt_Sigma.inv_X)
    beta <- crossprod(Xt_Sigma.inv_X.inv, crossprod(Sigma.inv_X, Y))
    
    # calc Xb
    fits <- X %*% beta
    # calc marginal residuals = (Y - Xb)
    residM <- as.vector(Y - fits)

    ### calculate PY
    PY <- crossprod(Sigma.inv, residM)
    
    # compute RSS
    YPY <- crossprod(Y, PY)
    RSS <- as.numeric(YPY/(n-k))
    # Sigma.inv_R <- crossprod(Sigma.inv, residM)
    # Rt_Sigma.inv_R <- crossprod(residM, Sigma.inv_R)
    # Rt_Sigma.inv_R <- crossprod(residM, PY)
    # RSS <- as.numeric(Rt_Sigma.inv_R/(n - k)) 

    # log likelihood
    logLik <- as.numeric(-0.5 * n * log(2 * pi * RSS) - sum(log(cholSigma.diag)) - 0.5 * YPY/RSS)
    # REML log likelihood; accounting for estimation of mean effects 
    logLikR <- as.numeric(logLik + 0.5 * k * log(2 * pi * RSS) - sum(log(diag(chol.Xt_Sigma.inv_X))))

    return(list(PY = PY, RSS = RSS, logLik = logLik, logLikR = logLikR, 
                Sigma.inv_X = Sigma.inv_X, Xt_Sigma.inv_X.inv = Xt_Sigma.inv_X.inv, 
                beta = as.numeric(beta), fits = fits, residM = residM))
    
}
