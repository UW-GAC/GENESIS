
## takes a null model and prepare specific arguments to streamline the testing
nullModelTestPrep <- function(nullmod){
    Y <- nullmod$fit$workingY
    X <- nullmod$model.matrix
    C <- nullmod$cholSigmaInv

    if (length(C) > 1) { ## n by n cholSigmaInv (may be Diagonal)
        if (is(C, "Matrix")) X <- Matrix(X, sparse = FALSE)
        CX <- crossprod(C, X)
        CXCXI <- tcrossprod(CX, chol2inv(chol(crossprod(CX))))
        # qrmod <- base::qr(CX)
        # Ytilde <- base::qr.resid(qrmod, as.matrix(crossprod(C, Y)))
        CY <- crossprod(C, Y)
        Ytilde <- CY - tcrossprod(CXCXI, crossprod(CY, CX))
        resid.PY <- C %*% Ytilde
        # resid <- tcrossprod(C, crossprod(nullmod$resid.marginal, C))

    } else { ## cholSigmaInv is a scalar
        CX <- C*X
        CXCXI <- tcrossprod(CX, chol2inv(chol(crossprod(CX))))
        # qrmod <- base::qr(CX)
        # Ytilde <- base::qr.resid(qrmod, as.matrix(C*Y))
        CY <- C*Y
        Ytilde <- CY - tcrossprod(CXCXI, crossprod(CY, CX))
        resid.PY <- C*Ytilde
        # resid <- nullmod$resid.marginal*C^2
    }

    # compute residual sum of squares under the null model
    RSS0 <- as.numeric(crossprod(Ytilde))

    #return(list(Ytilde = Ytilde, resid = resid, ))
    out <- list(resid.cholesky = Ytilde, resid.PY = resid.PY,
                prep_elements = list(CX = CX, CXCXI = CXCXI, RSS0 = RSS0))
    return(out)
}


##  adjust genotypes for correlation structure and fixed effects
calcGtilde <- function(nullmod, G){
    C <- nullmod$cholSigmaInv
    if(length(C) > 1){ # n by n cholSigmaInv (may be Diagonal)
        CG <- crossprod(C, G)
    }else{ # cholSigmaInv is a scalar
        CG <- C*G
    }

    # calculate Gtilde
    nrowG <- as.numeric(nrow(CG))
    ncolG <- as.numeric(ncol(CG))
    if(length(C) == 1 || nrowG*ncolG <= 2^31){
        Gtilde <- CG - tcrossprod(nullmod$CXCXI, crossprod(CG, nullmod$CX))
        # base::qr.resid(nullmod$qr, CG) # QR seems to be slower unexpectedly

    }else{
        # too large when G sparse; break into multiple blocks
        nblock <- ceiling(nrowG*ncolG/2^31)
        blocks <- unname(split(1:ncolG, cut(1:ncolG, nblock)))
        Gtilde <- list()
        for(i in 1:length(blocks)){
            Gtilde[[i]] <- as.matrix(CG[,blocks[[i]]] - tcrossprod(nullmod$CXCXI, crossprod(CG[,blocks[[i]]], nullmod$CX)))
        }
        Gtilde <- do.call(cbind, Gtilde)
    }

    return(Gtilde)
}

## adjust genotypes for correlation structure and fixed effects using fast approximation from SAIGE
## replace C = sqrt(Sigma^{-1}) with W^{1/2} (diagonal matrix)
calcGtildeFast <- function(nullmod, G, r = 1){    
    X <- nullmod$model.matrix
    W <- nullmod$W
    # W is the diagonal of a matrix
    WX <- W*X
    XWX.inv <- solve(crossprod(X,WX))
    # G - X(X'WX)^{-1}(X'WG) (formula from SAIGE)
    Gtilde <- G - tcrossprod(X, crossprod(crossprod(WX, G), XWX.inv))
    # multiply by r*sqrt(W) so that Gtilde'Gtilde = variance
    Gtilde <- r*sqrt(W)*Gtilde

    return(Gtilde)
}
