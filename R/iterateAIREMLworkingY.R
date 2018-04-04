
## fit here is an object return by the glm() function.
## eta is X %*% beta
.calcWorkingYnonGaussian <- function(y, eta, family){
    mu <- family$linkinv(eta) # exp(eta)/(1 + exp(eta)) for binomial
    # weights
    vmu <- family$variance(mu) # mu(1-mu) for binomial
    # inverse of g'(mu)
    gmuinv <- family$mu.eta(eta) # = vmu for canonical link
    # working vector
    Y <- eta + (y - mu)/gmuinv

    return(list(Y=Y, vmu=vmu, gmuinv=gmuinv))
}


.iterateAIREMLworkingY <- function(glm.mod, X, family, start = NULL, covMatList, AIREML.tol = 1e-6,
                                   drop.zeros = TRUE, max.iter = 100, verbose = TRUE){
    y <- glm.mod$y
    eta <- glm.mod$linear.predictors
    working.y <- .calcWorkingYnonGaussian(y, eta, family)
    newstart <- start
    Yreps <- 0
    
    repeat({
        Yreps <- Yreps + 1
        if(verbose) message("Computing Variance Component Estimates...")
        if(verbose) message(paste(paste("Sigma^2_",c(names(covMatList)),sep="", collapse="     "), "log-lik", "RSS", sep="     "))
        
        # estimate variance components
        vc.mod <- .runAIREMLother(Y=working.y$Y, X=X, start=newstart, covMatList=covMatList, 
                                  AIREML.tol=AIREML.tol, drop.zeros=drop.zeros, max.iter=max.iter, 
                                  verbose=verbose, vmu=working.y$vmu, gmuinv=working.y$gmuinv)
        
        if (vc.mod$allZero == TRUE) {
            message("All variance components estimated as zero, using glm...")
            break()
        }
        # update parameters
        if(verbose) message("Updating WorkingY Vector...")
        working.y <- .calcWorkingYnonGaussian(y, vc.mod$eta, family)
        
        # current variance component estimate
        newstart <- vc.mod$varComp
        newstart[vc.mod$zeroFLAG] <- AIREML.tol
        
        # test for convergence
        stat <- sqrt(sum((vc.mod$eta - eta)^2))
        if(verbose) message(paste("Checking for Convergence...", stat, sep = "\t"))
        eta <- vc.mod$eta
        if(stat < AIREML.tol){ break() }
        
        if(Yreps == max.iter){
            vc.mod$converged <- FALSE
            warning("Maximum number of iterations for workingY reached without convergence!")
            break()
        }
    })
    
    return(list(vc.mod = vc.mod, working.y = working.y))
    
}
