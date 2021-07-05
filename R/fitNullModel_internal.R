### A function that fits a null model, assuming complete data - will have to be wrapped using a different function that performs checks.
## or checks will be added... for a start - complete valid data. nrow(X) == length(y), dimensions of covMatList similarly appropriate. We
## do not deal with IDs, just indices, as everything is assumed to match.
## X is assumed to have an intercept.
## non-gaussian families can only be binomial and poisson.

## y - outcome vector
## X - data.frame or model.matrix
.fitNullModel <- function(y, X, covMatList = NULL, group.idx = NULL, family = gaussian(), start = NULL,
                          AIREML.tol = 1e-4, max.iter = 100, EM.iter = 0,
                          drop.zeros = TRUE, return.small = FALSE, verbose = TRUE){

    ### checks
    if(!is.null(covMatList)){
        covMatList = .setCovMatNames(covMatList)
        # if any Matrix objects; coerce all to Matrix objects (coerced ones are not sparse)
        covMatList <- .checkMatrixType(covMatList)
    }

    if (is.null(colnames(X))){
        colnames(X) <- paste0("X", 1:ncol(X))
    }

    ### Gaussian family
    if (family$family == "gaussian"){
        if (is.null(covMatList) & is.null(group.idx)) {
            # linear regression
            mod <- lm(y ~ -1 + X)
            out <- .nullModOutReg(y, X, mod, family)
        }
        if (is.null(covMatList) & !is.null(group.idx)){
            vc.mod <- .runWLSgaussian(y, X, group.idx = group.idx, start = start,
                                      AIREML.tol = AIREML.tol, max.iter = max.iter,
                                      EM.iter = EM.iter, verbose = verbose)
            out <- .nullModOutWLS(y, X, vc.mod = vc.mod, family = family, group.idx = group.idx)
        }
        if (!is.null(covMatList)){
            # LMM
            if (is.null(group.idx)) group.idx <- list(resid.var = 1:length(y))
            vc.mod <- .runAIREMLgaussian(y, X, start = start, covMatList = covMatList,
                                         group.idx = group.idx, AIREML.tol = AIREML.tol, drop.zeros = drop.zeros,
                                         max.iter = max.iter, EM.iter = EM.iter, verbose = verbose)
            out <- .nullModOutMM(y = y, workingY = y, X = X, vc.mod = vc.mod,
                                 family = family, covMatList = covMatList,
                                 group.idx = group.idx, drop.zeros = drop.zeros)
        }
    }

    ### Non-Gaussian family
    if (family$family != "gaussian"){
        # initial fit with glm
        mod <- glm(y ~ X, family = family)

        if (!is.null(covMatList)){ ## iterate between computing workingY and estimating VCs.
            iterate.out <- .iterateAIREMLworkingY(glm.mod = mod, X = X, family = family,
                                                  start = start, covMatList = covMatList, AIREML.tol = AIREML.tol,
                                                  drop.zeros = drop.zeros, max.iter = max.iter, EM.iter = EM.iter,
                                                  verbose = verbose)

            vc.mod <- iterate.out$vc.mod
            working.y <- iterate.out$working.y

      	    ## check whether all variance components were estimated as zero:
            if (vc.mod$allZero == TRUE){
                out <- .nullModOutReg(y, X, mod, family)
                out$zeroFLAG <- TRUE
            } else{
                out <- .nullModOutMM(y = y, workingY = working.y$Y, X = X, vc.mod = vc.mod,
                                     family = family, covMatList = covMatList,
                                     vmu = working.y$vmu, gmuinv = working.y$gmuinv, drop.zeros = drop.zeros)
            }
        } else{
            out <- .nullModOutReg(y, X, mod, family)
        }
    }

    out.class <- class(out)
    nullprep <- nullModelTestPrep(out)
    # Add to fit data frame.
    out$fit$resid.PY <- as.vector(nullprep$resid.PY)
    out$fit$resid.cholesky <- as.vector(nullprep$resid.cholesky)
    out <- c(out, nullprep$prep_elements)
    if (return.small) {
        out <- nullModelSmall(out)
    }
    class(out) <- out.class

    return(out)

}
