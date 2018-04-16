
## takes a null model and prepare specific arguments to streamline the testing
nullModelTestPrep <- function(nullmod){
    
    Y <- nullmod$workingY
    X <- Matrix(nullmod$model.matrix)
    
    if (nullmod$family$mixedmodel | (nullmod$family$family == "gaussian")){
        C <- nullmod$cholSigmaInv
        if (length(C) > 1) { ## n by n cholSigmaInv (may be Diagonal)
            CX <- crossprod(C, X)
            CXCXI <- tcrossprod(CX, chol2inv(chol(crossprod(CX))))
            CCXCXICX <- tcrossprod(tcrossprod(C, t(CXCXI)), CX)
            Ytilde <- crossprod(C, Y) - crossprod(CCXCXICX, Y)
            resid <- tcrossprod(C, t(Ytilde)) - tcrossprod(CCXCXICX, t(Ytilde))
        } else { ## cholSigmaInv is a scalar
            CX <- X * C 
            CXCXI <- tcrossprod(CX, chol2inv(chol(crossprod(CX))))
            CCXCXICX <- tcrossprod(CXCXI*C, CX)
            Ytilde <- C*Y - crossprod(CCXCXICX, Y)
            resid <- nullmod$resid.marginal/nullmod$varComp
        }
    }

    if (!nullmod$family$mixedmodel & (nullmod$family$family != "gaussian")){
        sigma <- sqrt(nullmod$varComp)
        C <- Diagonal(x=sigma)
        CX <- X * sigma
        CXCXI <- tcrossprod(CX, chol2inv(chol(crossprod(CX))))
        CCXCXICX <- tcrossprod(tcrossprod(C, t(CXCXI)), CX)
        Ytilde <- crossprod(C, Y) - crossprod(CCXCXICX, Y)
        resid <- nullmod$resid.marginal
    }

    return(list(Ytilde=Ytilde, resid=resid, CX=CX, CXCXI=CXCXI))
}


##  adjust genotypes for correlation structure and fixed effects
calcXtilde <- function(nullmod, G){
    
    if (nullmod$family$mixedmodel | (nullmod$family$family == "gaussian")){
        C <- nullmod$cholSigmaInv
        if (length(C) > 1) { ## n by n cholSigmaInv (may be Diagonal)
            M1 <- crossprod(C, G)
        } else { ## cholSigmaInv is a scalar
            M1 <- G * C
        }
    }

    if (!nullmod$family$mixedmodel & (nullmod$family$family != "gaussian")){
        sigma <- sqrt(nullmod$varComp)
        C <- Diagonal(x=sigma)
        M1 <- crossprod(C, G)
    }

    rm(G)
    Xtilde <- M1 - tcrossprod(nullmod$CXCXI, crossprod(M1, nullmod$CX))
    return(Xtilde)
}
