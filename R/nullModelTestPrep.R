
## takes a null model and prepare specific arguments to streamline the testing
nullModelTestPrep <- function(nullmod){ 
    Y <- nullmod$workingY
    X <- nullmod$model.matrix
    C <- nullmod$cholSigmaInv

    if (length(C) > 1) { ## n by n cholSigmaInv (may be Diagonal)
        if (is(C, "Matrix")) X <- Matrix(X)
        CX <- crossprod(C, X)
        CXCXI <- tcrossprod(CX, chol2inv(chol(crossprod(CX))))
        # qrmod <- base::qr(CX)
        # Ytilde <- base::qr.resid(qrmod, as.matrix(crossprod(C, Y)))
        CY <- crossprod(C, Y)
        Ytilde <- CY - tcrossprod(CXCXI, crossprod(CY, CX))
        resid <- C %*% Ytilde
        # resid <- tcrossprod(C, crossprod(nullmod$resid.marginal, C))
        
    } else { ## cholSigmaInv is a scalar
        CX <- C*X
        CXCXI <- tcrossprod(CX, chol2inv(chol(crossprod(CX))))
        # qrmod <- base::qr(CX)
        # Ytilde <- base::qr.resid(qrmod, as.matrix(C*Y))
        CY <- C*Y
        Ytilde <- CY - tcrossprod(CXCXI, crossprod(CY, CX))
        resid <- C*Ytilde
        # resid <- nullmod$resid.marginal*C^2
    }

    return(list(Ytilde = Ytilde, resid = resid, CX = CX, CXCXI = CXCXI))
    # return(list(Ytilde = Ytilde, resid = resid, CX = CX, CXCXI = CXCXI, qr = qrmod))
}


##  adjust genotypes for correlation structure and fixed effects
##  this replaces calcXtilde; changed the name to be less confusing; X is covariates and G is genotypes
calcGtilde <- function(nullmod, G){
    C <- nullmod$cholSigmaInv

    if(length(C) > 1){ # n by n cholSigmaInv (may be Diagonal)
        CG <- crossprod(C, G)
    }else{ # cholSigmaInv is a scalar
        CG <- C*G
    }

    # calculate Gtilde
    CG - tcrossprod(nullmod$CXCXI, crossprod(CG, nullmod$CX))
    # base::qr.resid(nullmod$qr, CG) # QR seems to be slower unexpectedly
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
