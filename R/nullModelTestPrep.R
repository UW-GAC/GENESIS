
## takes a null model and prepare specific arguments to streamline the testing
nullModelTestPrep <- function(nullmod){
    
    Y <- nullmod$workingY
    X <- nullmod$model.matrix
    
    if (nullmod$family$mixedmodel){  ## n by n cholSigmaInv
        C <- nullmod$cholSigmaInv
        CX <- crossprod(C, X)
    }

    if (!nullmod$family$mixedmodel & (nullmod$family$family != "gaussian")){
        sigma <- sqrt(nullmod$varComp)
        C <- Diagonal(x=sigma)
        CX <- X * sigma
    }

    if (!nullmod$family$mixedmodel & (nullmod$family$family == "gaussian")){
        if (nullmod$hetResid) {  ## cholSigmaInv is diagonal
            C <- diag(nullmod$cholSigmaInv)
        } else { ## family is "gaussian", cholSigmaInv is a scalar.
            C <- nullmod$cholSigmaInv
        }
        CX <- X * C      ## this is equal to crossprod(diag(C), X) when C is a vector  
    }

    CXCXI = tcrossprod(CX, chol2inv(chol(crossprod(CX))))
    
    if (nullmod$family$mixedmodel){
        CCXCXICX <- tcrossprod(tcrossprod(C, t(CXCXI)), CX)
        Ytilde <- crossprod(C, Y) - crossprod(CCXCXICX, Y)
        resid <- tcrossprod(C, t(Ytilde)) - tcrossprod(CCXCXICX, t(Ytilde))
    }
    
    if (!nullmod$family$mixedmodel & (nullmod$family$family != "gaussian")){
        CCXCXICX <- tcrossprod(tcrossprod(C, t(CXCXI)), CX)
        Ytilde <- crossprod(C, Y) - crossprod(CCXCXICX, Y)
        resid <- nullmod$resid.marginal
    }
    
    if (!nullmod$family$mixedmodel & (nullmod$family$family == "gaussian")){
        CCXCXICX <- tcrossprod(CXCXI*C, CX)
        Ytilde <- C*Y - crossprod(CCXCXICX, Y)
      if (nullmod$hetResid){
          resid <- as.vector(C*Ytilde - tcrossprod(CCXCXICX, t(Ytilde)))
      } else{
          resid <- nullmod$resid.marginal/nullmod$varComp
      }
    }
    
    return(list(Ytilde=Ytilde, resid=resid, CX=CX, CXCXI=CXCXI))
}


##  adjust genotypes for correlation structure and fixed effects
calcXtilde <- function(nullmod, G){
    
    if (nullmod$family$mixedmodel){  ## n by n cholSigmaInv
        M1 <- crossprod(nullmod$cholSigmaInv,G)
    }
    
    if (!nullmod$family$mixedmodel & (nullmod$family$family != "gaussian")){
        ## above we use cholSigmaInv, here we use Sigma, but the math is the same.
        ## why is this case different from the gaussian case below?
        C <- Diagonal(x=sqrt(nullmod$varComp))
        M1 <- crossprod(C,G)
    }

    if (!nullmod$family$mixedmodel & (nullmod$family$family == "gaussian")){  ## a diagonal or scalar cholSigmaInv
        if (length(nullmod$cholSigmaInv) > 1) { ## cholSigmaInv is diagonal
            C <- diag(nullmod$cholSigmaInv)
        } else { ## family is "gaussian", cholSigmaInv is a scalar.
            C <- nullmod$cholSigmaInv
        }
        M1 <- G*C
    }	
    
    rm(G)
    Xtilde <- M1 - tcrossprod(nullmod$CXCXI, crossprod(M1, nullmod$CX))
    return(Xtilde)
}
