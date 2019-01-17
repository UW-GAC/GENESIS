
### preparing output arguments for regression models that are not mixed.
.nullModOutReg <- function(y, X, mod, family, group.idx = NULL){
    family$mixedmodel <- FALSE
    
    if (family$family == "gaussian"){
        varComp <- summary(mod)$sigma^2
        cholSigmaInv <- sqrt(1/varComp)
        workingY <- drop(y)
    }  else{
        varComp <- NULL
        vmu <- family$variance(mod$fitted)
        cholSigmaInv <- Diagonal(x=sqrt(vmu))
        workingY <- .calcWorkingYnonGaussian(y, eta = mod$linear.predictors, family)$Y
    }
    
    varCompCov <- NULL
    hetResid <- FALSE 
    fixef <- as.data.frame(summary(mod)$coef)
    varNames <- colnames(X) 
    rownames(fixef) <- varNames
    betaCov <- vcov(mod, complete=FALSE)
    dimnames(betaCov) <- list(varNames, varNames)
    fitted.values <- mod$fitted.values
    resid.marginal <-  residuals(mod, type = "working")
    logLik <- as.numeric(logLik(mod))
    AIC <- AIC(mod)
    converged <- ifelse(family$family == "gaussian", TRUE, mod$converged)
    zeroFLAG <- NULL
    RSS <- ifelse(family$family == "gaussian", sum(resid.marginal^2)/varComp/(nrow(X) - ncol(X)), 1)
    
    out <- list(family = family, hetResid = hetResid, varComp = varComp,
                varCompCov = varCompCov, fixef = fixef, betaCov = betaCov, 
                fitted.values = fitted.values, resid.marginal = resid.marginal, 
                logLik = logLik, AIC = AIC, workingY = workingY, outcome = y, 
                model.matrix = X, group.idx = group.idx, cholSigmaInv = cholSigmaInv, 
                converged = converged, zeroFLAG = zeroFLAG, RSS = RSS)
    class(out) <- "GENESIS.nullModel"
    return(out)
}



.nullModOutWLS <- function(y, X, vc.mod, family, group.idx = NULL){
    family$mixedmodel <- FALSE
    
    if (is.null(names(group.idx))){
        group.names <- paste0("G", seq_along(group.idx))
    } else { # there are names
        group.names <- names(group.idx)
    }
    
    varComp <- vc.mod$varComp
    names(varComp) <- group.names
    varCompCov <- solve(vc.mod$AI)
    dimnames(varCompCov) <- list(group.names, group.names)
    
    hetResid <- TRUE
    varNames <- colnames(X) 
    ## cholSigmaInv <- Diagonal(x=sqrt(diag(vc.mod$Sigma.inv)))
    cholSigmaInv.diag <- sqrt(diag(vc.mod$Sigma.inv))
    cholSigmaInv <- Diagonal(x=cholSigmaInv.diag)
    
    RSS <- vc.mod$RSS
   
    ## betaCov <- as.matrix(RSS * chol2inv(chol(crossprod(crossprod(cholSigmaInv, X)))))
    betaCov <- as.matrix(RSS * chol2inv(chol(crossprod(cholSigmaInv.diag*X))))
    dimnames(betaCov) <- list(varNames, varNames)
    
    SE <- sqrt(diag(betaCov))
    Stat <- (vc.mod$beta/SE)^2
    pval <- pchisq(Stat, df = 1, lower.tail = FALSE)

    fixef <- data.frame(Est = vc.mod$beta, SE = SE, Stat = Stat, pval = pval)
    rownames(fixef) <- varNames
    
    fitted.values <- as.vector(vc.mod$fits)
    resid.marginal <-  vc.mod$residM
    logLik <- vc.mod$logLik
    logLikR <- vc.mod$logLikR
    AIC <- 2 * (ncol(X) + length(varComp)) - 2 * logLik
    
    workingY <- drop(y)
    
    resid.conditional <- workingY - drop(vc.mod$eta) ### should be the same as resid.marginal
    
    converged <- TRUE
    zeroFLAG <- NULL

    out <- list(family = family, hetResid = hetResid, varComp = varComp, varCompCov = varCompCov, 
                fixef = fixef, betaCov = betaCov, fitted.values = fitted.values, 
                resid.marginal = resid.marginal, resid.conditional = resid.conditional, 
                logLik = logLik, logLikR  = logLikR, AIC = AIC, workingY = workingY, 
                outcome = y, model.matrix = X, group.idx = group.idx, cholSigmaInv = cholSigmaInv, 
                converged = converged, zeroFLAG = zeroFLAG, RSS = RSS)
    class(out) <- "GENESIS.nullModel"
    return(out)
}



.nullModOutMM <- function(y, workingY, X, vc.mod, family, covMatList, group.idx = NULL, vmu = NULL, gmuinv = NULL, drop.zeros = TRUE){
    n <- nrow(X)
    m <- length(covMatList)
    
    if (!is.null(group.idx)){
        g <- length(group.idx)
        if (is.null(names(group.idx))){
            group.names <- paste0("G", 1:g)
        } else{
            group.names <- names(group.idx)
        }
    }
    if (is.null(group.idx)){
        if (family$family == "gaussian"){
            group.idx <- list(E = 1:n)
            g <- 1
            group.names <- "E" 
        } else{
            g <- 0
            group.names <- NULL
        }
    }
    
    if(is.null(names(covMatList))){
        names(covMatList) <- paste0("M",1:m)
    }

    matrix.names <- names(covMatList)
    
    family$mixedmodel <- TRUE
    varComp <- vc.mod$varComp
    hetResid <- (length(group.idx) > 1)

    
    names(varComp) <- paste("V_",c(names(covMatList),group.names),sep="")
    varCompCov <- matrix(NA, nrow=(m+g), ncol=(m+g))
    colnames(varCompCov) <- paste("V_",c(names(covMatList),group.names),sep="")
    rownames(varCompCov) <- paste("V_",c(names(covMatList),group.names),sep="")

    
    if(drop.zeros){
        varCompCov[!vc.mod$zeroFLAG, !vc.mod$zeroFLAG] <- solve(vc.mod$AI)
    }else{
        varCompCov <- solve(vc.mod$AI)
    }
    
    sq <- .computeSigmaQuantities(varComp = varComp, covMatList = covMatList, group.idx = group.idx, vmu = vmu, gmuinv = gmuinv)
    # Cholesky Decomposition
    cholSigmaInv <- t(chol(sq$Sigma.inv))
    dimnames(cholSigmaInv) <- list(colnames(covMatList[[1]]), colnames(covMatList[[1]]))
    

    varNames <- colnames(X) 
    
    RSS <- ifelse(family$family == "gaussian", vc.mod$RSS, 1)

    X.tmp <- if (is(cholSigmaInv, "Matrix")) Matrix(X) else X
    betaCov <- as.matrix(RSS * chol2inv(chol(crossprod(crossprod(cholSigmaInv, X.tmp)))))
    rm(X.tmp)
    #betaCov <- as.matrix(RSS * chol2inv(chol(crossprod(crossprod(cholSigmaInv, X)))))

    dimnames(betaCov) <- list(varNames, varNames)
    
    SE <- sqrt(diag(betaCov))
    Stat <- (vc.mod$beta/SE)^2
    pval <- pchisq(Stat, df = 1, lower.tail = FALSE)

    fixef <- data.frame(Est = vc.mod$beta, SE = SE, Stat = Stat, pval = pval)
    rownames(fixef) <- varNames
    
    fitted.values <- as.vector(vc.mod$fits)
    resid.marginal <-  vc.mod$residM
    logLik <- vc.mod$logLik
    logLikR <- vc.mod$logLikR
    AIC <- 2 * (ncol(X) + length(varComp)) - 2 * logLik
    
    workingY <- drop(workingY)
    
    resid.conditional <- workingY - drop(vc.mod$eta)
    
    converged <- vc.mod$converged
    zeroFLAG <- vc.mod$zeroFLAG

    out <- list(family = family, hetResid = hetResid, varComp = varComp, 
                varCompCov = varCompCov, fixef = fixef, betaCov = betaCov, 
                fitted.values = fitted.values, resid.marginal =resid.marginal, 
                resid.conditional = resid.conditional, logLik = logLik, 
                logLikR = logLikR, AIC = AIC, workingY = workingY, outcome = y,
                model.matrix = X, group.idx = group.idx, cholSigmaInv = cholSigmaInv, 
                converged = converged, zeroFLAG = zeroFLAG, RSS = RSS)
    class(out) <- "GENESIS.nullMixedModel"
    return(out)
}
