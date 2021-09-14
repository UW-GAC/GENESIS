
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
    # diagV
    diagV <- as.vector(vmu)/as.vector(gmuinv)^2

    return(list(Y=Y, diagV=diagV))
}


.iterateAIREMLworkingY <- function(glm.mod, X, family, start = NULL, covMatList, AIREML.tol = 1e-4,
                                   drop.zeros = TRUE, max.iter = 100, EM.iter = 0, verbose = TRUE){
    y <- glm.mod$y
    eta <- glm.mod$linear.predictors
    working.y <- .calcWorkingYnonGaussian(y, eta, family)
    newstart <- start
    Yreps <- 0

    repeat({
        Yreps <- Yreps + 1

        # estimate variance components
        vc.mod <- .runAIREMLother(Y = working.y$Y, X = X, start = newstart, covMatList = covMatList,
                                  diagV = working.y$diagV, AIREML.tol = AIREML.tol, drop.zeros = drop.zeros,
                                  max.iter = max.iter, EM.iter = EM.iter, verbose = verbose)

        if (vc.mod$allZero == TRUE) {
            message("All variance components estimated as zero, using glm...")
            break()
        }

        ### check for convergence
        if(sqrt(sum((vc.mod$eta - eta)^2)) < AIREML.tol){
            converged <- TRUE
            (break)()
        }else{
            # check if exceeded the number of iterations
            if(Yreps == max.iter){
                vc.mod$converged <- FALSE
                warning("Maximum number of iterations for workingY reached without convergence!")
                (break)()
            }else{
                # update starting values to current variance component estimates
                newstart <- vc.mod$varComp
                newstart[vc.mod$zeroFLAG] <- AIREML.tol
                # update workingY
                if(verbose) message("Updating WorkingY Vector...")
                working.y <- .calcWorkingYnonGaussian(y, vc.mod$eta, family)
                # update eta
                eta <- vc.mod$eta
            }
        }
    })

    return(list(vc.mod = vc.mod, working.y = working.y))

}
