### A function that fits a null model, assuming complete data - will have to be wrapped using a different function that performs checks. 
## or checks will be added... for a start - complete valid data. nrow(X) == length(y), dimensions of covMatList similarly appropriate. We
## do not deal with IDs, just indices, as everything is assumed to match. 
## X is assumed to have an intercept. 
## non-gaussian families can only be binomial and poisson. 

## y - outcome vector
## X - data.frame or model.matrix
.fitNullModel <- function(y, X, covMatList = NULL, group.idx = NULL, family = "gaussian", start = NULL,
                          AIREML.tol = 1e-6, max.iter = 100, drop.zeros = TRUE, verbose = TRUE){
    
    if(!is.null(covMatList)){
        if (!is.list(covMatList)){
            covMatList <- list(A = covMatList)
        }
        # coerce to Matrix objects. should get "dspMatrix" (packed symmetric matrix)
        covMatList <- lapply(covMatList, function(x) {
            x <- Matrix(x)
            if (is(x, "symmetricMatrix") & !is(x, "sparseMatrix")) x <- pack(x)
            return(x)
        })
    }

    if (is.null(colnames(X))){
        colnames(X) <- paste0("X", 1:ncol(X))
    }
    
    if(is.character(family)){
        family <- get(family)
    }
    if(is.function(family)){
        family <- family()
    }
    if(is.null(family$family)){
        stop("'family' not recognized")
    }
    if (!is.element(family$family, c("gaussian", "binomial", "poisson"))){
        stop("family must be one of gaussian, binomial, or poisson")
    }

    ## save original model.matrix for output, in case we need to convert to Matrix
    X.mm <- X
    
    if (family$family == "gaussian"){
        if (is.null(covMatList) & is.null(group.idx)) {
            mod <- lm(y ~ -1 + X)  ## prepare output based on that. 
            out <- .nullModOutReg(y, X, mod, family)
        }
        if (is.null(covMatList) & !is.null(group.idx)){
            X <- Matrix(X)
            vc.mod <- .runWLSgaussian(y, X, group.idx = group.idx, start = start, 
                                      AIREML.tol = AIREML.tol, max.iter = max.iter, verbose = verbose)
            out <- .nullModOutWLS(y, X, vc.mod = vc.mod, family = family, group.idx = group.idx)
        }
        if (!is.null(covMatList)){
            X <- Matrix(X)
            if (is.null(group.idx)) group.idx <- list(resid.var = 1:length(y))
            vc.mod <- .runAIREMLgaussian(y, X, start = start, covMatList = covMatList, 
                                         group.idx = group.idx, AIREML.tol = AIREML.tol, drop.zeros = drop.zeros,  
                                         max.iter = max.iter, verbose = verbose)
            out <- .nullModOutMM(y = y, workingY = y, X = X, vc.mod = vc.mod, 
                                 family = family, covMatList = covMatList, 
                                 group.idx = group.idx, drop.zeros = drop.zeros)
        }
    } 
    if (family$family != "gaussian"){ # separate condition instead of "else" for readability. 
        mod <- glm(y ~ X, family = family)
        
        if (!is.null(covMatList)){ ## iterate between computing workingY and estimating VCs. 
            X <- Matrix(X)
            iterate.out <- .iterateAIREMLworkingY(glm.mod = mod, X = X, family = family, 
                                                  start = start, covMatList = covMatList, AIREML.tol = AIREML.tol,
                                                  drop.zeros = drop.zeros, max.iter = max.iter, verbose = verbose)
            
            vc.mod <- iterate.out$vc.mod
            working.y <- iterate.out$working.y
            
      	    ## check whether all variance components were estimated as zero:
            if (vc.mod$allZero == TRUE){
                out <- .nullModOutReg(y, X, mod, family)
                out$zeroFLAG <- TRUE
            } else{
                out <- .nullModOutMM(y = y, workingY = working.y$Y, X = X, 
                                     vc.mod = vc.mod, family = family, covMatList = covMatList, 
                                     vmu = working.y$vmu, gmuinv = working.y$gmuinv, drop.zeros = drop.zeros)
            }	
        } else{
            out <- .nullModOutReg(y, X, mod, family)
        }
    }

    out.class <- class(out)
    out$model.matrix <- Matrix(out$model.matrix)
    nullprep <- nullModelTestPrep(out)
    out <- c(out, nullprep)
    out$model.matrix <- X.mm # preserve original model.matrix object in output
    class(out) <- out.class
    
    return(out)
    
}
