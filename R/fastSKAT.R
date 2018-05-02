
nullModelFastSKATTestPrep <- function(nullmod, threshold = 1e-10){
    
    Y <- nullmod$workingY
    X <- nullmod$model.matrix
    C <- nullmod$cholSigmaInv

    ## C can also be a scalar (simple linear regression) - adapt code  
    if (is.null(dim(C)) & length(C) == 1){
        cholSigma <- 1/C
        SigmaInv <- C^2
        qr <- qr(X/cholSigma)
    } else{
        
        SIGMA <- chol2inv(t(C))
        SIGMA[abs(SIGMA) < threshold] <- 0
        SIGMA <- Matrix(SIGMA)
        cholSigma <- t(chol(SIGMA))
        SigmaInv <- tcrossprod(C)

        X <- Matrix(X)

        qr <- qr(as.matrix(solve(cholSigma,X)))
    }

    out <- list(cholSigma = cholSigma, SigmaInv = SigmaInv, qr = qr, resid = nullmod$resid.marginal)
    return(out)
}


nullModelFastSKAT <- function(null.model, threshold = 1e-10){
    nm.class <- class(null.model)
    out <- nullModelFastSKATTestPrep(null.model, threshold=threshold)
    null.model[c("cholSigmaInv", "Ytilde", "resid", "CX", "CXCXI")] <- NULL
    null.model <- c(null.model, out)
    class(null.model) <- c(nm.class, "GENESIS.nullModelFastSKAT")
    return(null.model)
}



.testVariantSetFastSKAT <- function(nullmod, G, weights, rho=0,
                                    method=c("ssvd","lanczos","satterthwaite"),
                                    neig=100, tr2.sample.size=500, q=NULL,
                                    convolution.method=c("saddlepoint","integration"),
                                    remainder.underflow=c("warn","missing","error")){
    
    if (!requireNamespace("bigQF")) stop("package 'bigQF' must be installed to run fastSKAT")
    if (!is(nullmod, "GENESIS.nullModelFastSKAT")) stop("run 'nullModelFastSKAT' to prepare the null model for this test")

    mindim <- min(dim(G))
    if (mindim < neig) neig <- mindim

    method <- match.arg(method)
    convolution.method <- match.arg(convolution.method)
    remainder.underflow <- match.arg(remainder.underflow)
    
    if (rho == 0){
        sparseG <- Matrix(G, sparse = TRUE) %*% Diagonal(x = weights) 
    } else{
        rho.mat <- matrix(rho, nrow = length(weights), ncol = length(weights))
        diag(rho.mat) <- 1
        cholRhoMat <- t(chol(rho.mat, pivot=TRUE))
        sparseG <- Matrix(G, sparse = TRUE) %*% Diagonal(x = weights)  %*% cholRhoMat
    }
    
    if (is.null(dim(nullmod$cholSigma)) & length(nullmod$cholSigma) == 1){
        rval <- list(mult = function(X){
            base::qr.resid(nullmod$qr, as.matrix((sparseG %*% X)/nullmod$cholSigma))	
        }, tmult = function(X){
            crossprod(sparseG, base::qr.resid(nullmod$qr,X)/nullmod$cholSigma)
        }, 
        trace = NULL,
        ncol = ncol(G),
        nrow = nrow(G),
        Q = function(){
            stdres <- nullmod$SigmaInv * nullmod$resid
            s = crossprod(sparseG, stdres)
            sum(s^2)
        })
        
        class(rval) <- c( "matrixfree") 
	
    } else{
        rval <- list(mult = function(X){
            base::qr.resid(nullmod$qr, as.matrix(solve(nullmod$cholSigma, (sparseG %*% X))))	
        }, tmult = function(X){
            crossprod(sparseG, solve(t(nullmod$cholSigma), base::qr.resid(nullmod$qr,X)))
        }, 
        trace = NULL,
        ncol = ncol(G),
        nrow = nrow(G),
        Q = function(){
            stdres <- nullmod$SigmaInv %*% nullmod$resid
            s = crossprod(sparseG, stdres)
            sum(s^2)
        })

        class(rval) <- c("famSKAT_genesis", "famSKAT", "matrixfree") 
        
    }
    
    pval <- bigQF::pQF(rval$Q(), rval,  method = method, neig = neig, tr2.sample.size = tr2.sample.size, q = q, convolution.method = convolution.method, remainder.underflow = remainder.underflow) ### what object should be the second entry of pQF??

    out <- list(pval=pval)
    return(out)
}
