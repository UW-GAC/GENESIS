
## takes a null model and prepare specific arguments to streamline the testing
nullModelTestPrep <- function(nullmod){
    
    Y <- nullmod$workingY
    X <- nullmod$model.matrix
    
    C <- nullmod$cholSigmaInv
    if (length(C) > 1) { ## n by n cholSigmaInv (may be Diagonal)
        if (is(C, "Matrix")) X <- Matrix(X)

        qrmod <- base::qr(crossprod(C, X))
        Ytilde <- base::qr.resid(qrmod, as.matrix(crossprod(C, Y)))
        resid <- C %*% Ytilde

        # CX <- crossprod(C, X)
        # CXCXI <- tcrossprod(CX, chol2inv(chol(crossprod(CX))))
###        CCXCXICX <- tcrossprod(tcrossprod(C, t(CXCXI)), CX)
###        Ytilde <- crossprod(C, Y) - crossprod(CCXCXICX, Y)
###        resid <- tcrossprod(C, t(Ytilde)) - tcrossprod(CCXCXICX, t(Ytilde))
        # CY <- crossprod(C, Y)
        # Ytilde <- CY - tcrossprod(CXCXI, crossprod(CY, CX))
        # resid <- tcrossprod(C, crossprod(nullmod$resid.marginal, C))
        
    } else { ## cholSigmaInv is a scalar

        qrmod <- base::qr(C*X)
        Ytilde <- base::qr.resid(qrmod, as.matrix(C*Y))
        resid <- C*Ytilde

        # CX <- X * C
        # CXCXI <- tcrossprod(CX, chol2inv(chol(crossprod(CX))))
###        CCXCXICX <- tcrossprod(CXCXI*C, CX)
###        Ytilde <- C*Y - crossprod(CCXCXICX, Y)
        # CY <- C*Y
        # Ytilde <- CY - tcrossprod(CXCXI, crossprod(CY, CX))
        # resid <- nullmod$resid.marginal*C^2
    }
    return(list(Ytilde = Ytilde, resid = resid, qr = qrmod))
    # return(list(Ytilde=Ytilde, resid=resid, CX=CX, CXCXI=CXCXI))
}


##  adjust genotypes for correlation structure and fixed effects
##  this replaces calcXtilde; changed the name to be less confusing; X is covariates and G is genotypes
calcGtilde <- function(nullmod, G){
    C <- nullmod$cholSigmaInv
    if(length(C) > 1){ # n by n cholSigmaInv (may be Diagonal)
        base::qr.resid(nullmod$qr, as.matrix(crossprod(C, G)))
    }else{
        base::qr.resid(nullmod$qr, as.matrix(C*G))
    }
}

##  adjust genotypes for correlation structure and fixed effects
# calcXtilde <- function(nullmod, G){
#     C <- nullmod$cholSigmaInv
#     if (length(C) > 1) { ## n by n cholSigmaInv (may be Diagonal)
#         M1 <- crossprod(C, G)
#     } else { ## cholSigmaInv is a scalar
#         M1 <- G * C
#     }

#     rm(G)
#     Xtilde <- M1 - tcrossprod(nullmod$CXCXI, crossprod(M1, nullmod$CX))
#     return(Xtilde)
# }
