
### add heterogeneous variants - covMatList == NULL
# n - number of people 
# g - the number of groups (1 to more). if g == 1 group.idx is ignored. 
.runWLSgaussian <- function (Y, X, group.idx, start, AIREML.tol,
                             max.iter, verbose){
    n <- length(Y)
    g <- length(group.idx)
    if (g <= 1) stop("group.idx must have length > 1")
    ## initializing parameters
    sigma2.p <- var(Y)
    AIREML.tol <- AIREML.tol * sigma2.p
    val <- 2 * AIREML.tol
    if (is.null(start)) {
        sigma2.k <- rep(sigma2.p,  g)
    } else {
        sigma2.k <- as.vector(start)
    }
    reps <- 0
    
    ## starting iteration
    repeat ({
        reps <- reps + 1
        
        if (reps > max.iter) {
            converged <- FALSE
            warning("Maximum number of iterations reached without convergence!")
            (break)()
        }

        
        ### Sigma is the matrix of variances. All non-diagonals are zero. 
        ## initialize Sigma (we don't need Sigma = only its diagonal!)
        diagSigma <- rep(0, n)
        
        ## set the values of the diagonal to be the group-specific variances. 
        for (i in 1:g) {
            diagSigma[group.idx[[i]]] <- sigma2.k[i]
        }
        
        ## just the diagonal - squared root of diagonals, and inverse of diagonals of Sigma.
        cholSigma.diag <- sqrt(diagSigma)
        Sigma.inv <- Diagonal(x=1/diagSigma)
        
        lq <- .calcLikelihoodQuantities(Y, X, Sigma.inv, cholSigma.diag, PPY=TRUE)
        
        if (verbose)  print(c(sigma2.k, lq$logLikR, lq$RSS))

        
        ## Updating variances and calculating their covariance matrix
        if (reps > 1) {
            
            score.AI <- .calcAIhetvars(lq$PY, lq$PPY, group.idx,
                                       Sigma.inv = Sigma.inv, Sigma.inv_X = lq$Sigma.inv_X, Xt_Sigma.inv_X.inv = lq$Xt_Sigma.inv_X.inv)
            score <- score.AI$score
            AI <- score.AI$AI
            AIinvScore <- solve(AI, score)

            sigma2.kplus1 <- sigma2.k + AIinvScore
            
            tau <- 1
            while (!all(sigma2.kplus1 >= 0)) {
                tau <- 0.5 * tau
                sigma2.kplus1 <- sigma2.k + tau * AIinvScore

            }
            val <- sqrt(sum((sigma2.kplus1 - sigma2.k)^2))
            sigma2.k <- sigma2.kplus1

            if (val < AIREML.tol) {
            	converged <- TRUE
                (break)()
            }
        }
        else { # first rep
            sigma2.kplus1 <- rep(NA, g)

            for (i in 1:g) {
                ### sigma2.kplus1[i] <- (1/n) * (sigma2.k[i]^2 * crossprod(lq$PY[group.idx[[i]]]) + n *sigma2.k[i] - sigma2.k[i]^2 * sum(diag(lq$P)[group.idx[[i]]]))
                covMati <- Diagonal( x=as.numeric( 1:n %in% group.idx[[i]] ) )
                trPi.part1 <- sum(diag(Sigma.inv)[ group.idx[[i]] ] )
                trPi.part2 <- sum(diag( (crossprod( lq$Sigma.inv_X, covMati) %*% lq$Sigma.inv_X) %*% lq$Xt_Sigma.inv_X.inv ))
                trPi <- trPi.part1 - trPi.part2
                sigma2.kplus1[i] <- as.numeric((1/n)*(sigma2.k[i]^2*crossprod(lq$PY[group.idx[[i]]]) + n*sigma2.k[i] - sigma2.k[i]^2*trPi ))
            }

            sigma2.k <- sigma2.kplus1
        }
    })
    
    ## after convergence, updated sigma again
    for (i in 1:g) {
        diagSigma[group.idx[[i]]] <- sigma2.k[i]
    }
    
    ## just the diagonal - squared root of diagonals, and inverse of diagonals of Sigma.
    cholSigma.diag <- sqrt(diagSigma)
    Sigma.inv <- Diagonal(x=1/diagSigma)
    
    lq <- .calcLikelihoodQuantities(Y, X, Sigma.inv, cholSigma.diag, PPY=TRUE)
    score.AI <- .calcAIhetvars(lq$PY, lq$PPY, group.idx,
                               Sigma.inv = Sigma.inv, Sigma.inv_X = lq$Sigma.inv_X, Xt_Sigma.inv_X.inv = lq$Xt_Sigma.inv_X.inv)
    AI <- score.AI$AI


    eta <- as.numeric(lq$fits)
    return(list(varComp = sigma2.k, AI = AI, converged = converged,
                Sigma.inv = Sigma.inv, beta = lq$beta, 
                residM = lq$residM, fits = lq$fits, eta = eta, logLikR = lq$logLikR,
                logLik = lq$logLik, RSS = lq$RSS))
}



