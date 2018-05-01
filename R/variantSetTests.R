
## function that gets an n \times p matrix of p genotypes of n individuals, and a null model, and tests the genotypes associations with the outcomes. 
## Genetic data are always assumed complete. 
## Types of tests: 
## Variant set: SKAT, burden, SKAT-O. Multiple types of p-values. Default: Davies with Kuonen if does not converge. 


testVariantSet <- function(nullmod, G, weights, test = c("Burden", "SKAT", "SMMAT", "fastSKAT"), 
                           burden.test = c("Score", "Wald"),  rho = 0,
                           pval.method = c("davies", "kuonen", "liu"), 
                           return.scores = FALSE, return.scores.cov = FALSE,
                           ...){

    test <- match.arg(test)
    burden.test <- match.arg(burden.test)
    pval.method <- match.arg(pval.method)

    if (is.matrix(G)) G <- Matrix(G)

    if (test == "SKAT") {
        out <- .testVariantSetSKAT(nullmod, G, weights, rho, pval.method, 
                                   return.scores, return.scores.cov)
    }    
    if (test == "Burden") {
        out <- .testVariantSetBurden(nullmod, G, weights, burden.test)
    }
    if (test == "SMMAT") {
        out <- .testVariantSetSMMAT(nullmod, G, weights, pval.method)
    }
    if (test == "fastSKAT") {
        out <- .testVariantSetFastSKAT(nullmod, G, weights, rho, ...)
    }
    return(out)
}



## create the burden score, than calls the appropriate single variant test function. 
## can easily implement GxE interaction with the burden score... later!
.testVariantSetBurden <- function(nullmod, G, weights, burden.test){
    
    burden <- Matrix(colSums(t(G) * weights))
    
    Xtilde <- calcXtilde(nullmod, burden)
    
    if (burden.test == "Score") {
        out <- .testGenoSingleVarScore(Xtilde, G = burden, resid = nullmod$resid) 
    }
    if (burden.test == "Wald"){
        out <- .testGenoSingleVarWald(Xtilde, Ytilde = nullmod$Ytilde,
                                      n = length(nullmod$Ytilde), k = ncol(nullmod$model.matrix))
    }
    return(out)
}



.testVariantSetSKAT <- function(nullmod, G, weights, rho = 0, pval.method, 
                                return.scores = FALSE, return.scores.cov = FALSE){

    U <- as.vector(crossprod(G, nullmod$resid))
    G <- calcXtilde(nullmod, G)
    if (length(rho) == 1) {
        out <- .runSKATTest(scores = U, geno.adj = G,
                            weights = weights, rho = rho, pval.method = pval.method,
                            optimal = FALSE)
    } else {
        ## SKAT-O
        out <- .runSKATTest(scores = U, geno.adj = G,
                            weights = weights, rho = rho, pval.method = pval.method,
                            optimal = TRUE)
    }
    return(out)
}



.testVariantSetSMMAT <- function(nullmod, G, weights, pval.method) {
    G <- t(t(G) * weights)
    U <- as.vector(crossprod(G, nullmod$resid))
    G <- calcXtilde(nullmod, G)
    V <- crossprod(G)
    burden.scores <- sum(U)
    burden.distMat <- sum(V)
    burden.pval <- pchisq(burden.scores^2/burden.distMat, df=1, lower.tail=FALSE)
    V.rowSums <- rowSums(V)
    U <- U - V.rowSums * burden.scores / burden.distMat
    V <- V - tcrossprod(V.rowSums) / burden.distMat
    if(mean(abs(V)) < sqrt(.Machine$double.eps)) return(list(pval_burden=burden.pval, pval_SMMAT=burden.pval, err=0))
    Q <- sum(U^2)
    # lambda for p value calculation
    lambda <- eigen(V, only.values = TRUE, symmetric=TRUE)$values
    lambda <- lambda[lambda > 0]
    pv <- .calcPval(Q, lambda, pval.method)
    pval <- pv["pval"]
    err <- pv["err"]
    out.pval <- tryCatch(pchisq(-2*log(burden.pval)-2*log(pval), df=4, lower.tail = FALSE), error = function(e) { NA })
    if(is.na(out.pval)) {
        err <- 1
        out.pval <- burden.pval
    }
    out <- list(pval_burden=burden.pval, pval_SMMAT=out.pval, err=err)
    return(out)
}



.runSKATTest <- function(scores, geno.adj, weights, rho, pval.method, optimal){
    # covariance of scores
    V <- crossprod(geno.adj)
    
    # vectors to hold output
    nrho <- length(rho)
    out.Q <- rep(NA, nrho); names(out.Q) <- paste("Q", rho, sep="_")
    out.pval <- rep(NA, nrho); names(out.pval) <- paste("pval", rho, sep="_")
    out.err <- rep(NA, nrho); names(out.err) <- paste("err", rho, sep="_")

    # for SKAT-O need to save lambdas
    if(optimal){ lambdas <- vector("list", nrho) }

    for(i in 1:nrho){
        if(rho[i] == 0){
            # Variance Component Test
            Q <- sum((weights*scores)^2) # sum[(w*scores)^2]  # for some reason SKAT_emmaX divides this by 2
            distMat <- weights*t(weights*V)  # (weights) V (weights) = (weights) X' P X (weights)

        }else if(rho[i] == 1){
            # Burden Test
            Q <- sum(weights*scores)^2 # (sum[w*scores])^2  # for some reason SKAT_emmaX divides this by 2
            distMat <- crossprod(weights,crossprod(V, weights)) # weights^T V weights

        }else if(rho[i] > 0 & rho[i] < 1){
            rhoMat <- matrix(rho[i], nrow=length(scores), ncol=length(scores)); diag(rhoMat) <- 1
            cholRhoMat <- t(chol(rhoMat, pivot=TRUE))
            Q <- crossprod(crossprod(weights*cholRhoMat,scores)) # scores' (weights) (rhoMat) (weights) scores
            distMat <- crossprod(cholRhoMat, crossprod(weights*t(weights*V), cholRhoMat)) # (cholRhoMat) (weights) X' P X (weights) (cholRhoMat)
        }

        # lambda for p value calculation
        lambda <- eigen(distMat, only.values = TRUE, symmetric=TRUE)$values
        lambda <- lambda[lambda > 0]
        if(optimal){ lambdas[[i]] <- lambda }

        # p value calculation
        if(length(scores) == 1){
            pval <- pchisq(as.numeric(Q/distMat), df=1, lower.tail=FALSE)
            err <- 0

        }else{
            pv <- .calcPval(Q, lambda, pval.method)
            pval <- pv["pval"]
            err <- pv["err"]
        }

        # update results
        out.Q[i] <- Q
        out.pval[i] <- pval
        out.err[i] <- err
    }
    out <- as.list(c(out.Q, out.pval, out.err))

    # SKAT-O
    if(optimal){
        if(length(scores) == 1){
            # pvalue is the same for all rhos
            out2 <- list(min.pval=out.pval[1], opt.rho=NA, pval_SKATO=out.pval[1])

        }else{
            # find the minimum p-value
            minp <- min(out.pval)
            out2 <- list(min.pval=minp, opt.rho=rho[which.min(out.pval)])

            # get qmin(rho); i.e. the (1-minp)th percentile of dist of each Q
            qmin <- rep(NA, nrho)
            for(i in 1:nrho){
                qmin[i] <- skatO_qchisqsum(minp, lambdas[[i]])
            }

            # calculate other terms
            Z <- t(t(geno.adj)*weights)
            zbar <- rowMeans(Z)
            zbarTzbar <- sum(zbar^2)
            M <- tcrossprod(zbar)/zbarTzbar
            ZtImMZ <- crossprod(Z, crossprod(diag(nrow(M)) - M, Z))
            lambda.k <- eigen(ZtImMZ, symmetric = TRUE, only.values = TRUE)
            lambda.k <- lambda.k$values[lambda.k$values > 0]
            mua <- sum(lambda.k)
            sum.lambda.sq <- sum(lambda.k^2)
            sig2a <- 2*sum.lambda.sq
            trMatrix <- crossprod(crossprod(Z,crossprod(M,Z)),ZtImMZ)
            sig2xi <- 4*sum(diag(trMatrix))
            kera <- sum(lambda.k^4)/sum.lambda.sq^2 * 12
            ldf <- 12/kera

            # calculate tau(rho)
            tau <- ncol(Z)^2*rho*zbarTzbar + (1-rho)*sum(crossprod(zbar, Z)^2)/zbarTzbar

            # find min{(qmin(rho)-rho*chisq_1)/(1-rho)} with integration							
            otherParams <- c(mu = mua, degf = ldf, varia = sig2a+sig2xi)
            # integrate
            re <- tryCatch({
                integrate(integrateFxn, lower = 0, upper = 40, subdivisions = 2000, qmin = qmin, otherParams = otherParams, tau = tau, rho = rho, abs.tol = 10^-25)
            }, error=function(e) NA)
            out2[["pval_SKATO"]] <- 1-re[[1]]
        }
        
        # update results
        out <- c(out, out2)
    }

    # return results
    return(out)	
}


.calcPval <- function(Q, lambda, pval.method) {
    if(!requireNamespace("survey")) stop("package 'survey' must be installed to calculate p-values for SKAT")
    if(!requireNamespace("CompQuadForm")) stop("package 'CompQuadForm' must be installed to calculate p-values for SKAT")
    
    err <- 0
    if(pval.method == "kuonen"){
        pval <- survey:::saddle(x = Q, lambda = lambda)
        err <- ifelse(is.na(pval), 1, 0)

    }else if(pval.method == "davies"){
        tmp <- suppressWarnings(CompQuadForm::davies(q = Q, lambda = lambda, acc = 1e-06))
        pval <- tmp$Qq
        if((tmp$ifault > 0) | (pval <= 0) | (pval >= 1)) {
            pval <- survey:::saddle(x = Q, lambda = lambda)
        }
        err <- ifelse(is.na(pval), 1, 0)

    }else if(pval.method == "liu"){
        pval <- CompQuadForm::liu(q = Q, lambda = lambda)
        err <- 0
    }
    
    if(err > 0){
        pval <- CompQuadForm::liu(q = Q, lambda = lambda)
    }

    return(c(pval=pval, err=err))
}


# function to calculate q_min value
# basically a qchisqsum() function that takes the quantile/percentile and the lambda values
# matches the first 2 moments and the kurtosis
# based upon liu et al (2009) paper
skatO_qchisqsum <- function(p, lambdas){
    mu <- sum(lambdas)
    sum.lambda.sq <- sum(lambdas^2)

    s1 <- sum(lambdas^3)/(sum.lambda.sq^(3/2))
    s2 <- sum(lambdas^4)/(sum.lambda.sq^2)
    if(s1^2 > s2){
    	a <- 1/(s1-sqrt(s1^2-s2))
    	d <- s1*a^3 - a^2
    	l <- a^2 - 2*d
    }else{ # s1^2 <= s2
        l <- 1/s2 # in liu et al, this is l=1/s1^2; matches kurtosis instead of skewness to improve tail prob estimates
    }  	
    
    qmin <- qchisq(1-p, df=l)
    pval <- (qmin - l)/sqrt(2*l) * sqrt(2*sum.lambda.sq) + mu
    
    return(pval)
}


## function to integrate; the first term of the optimal integrand
# it's a non-central sum of weighted chi-squares
integrateFxn <- function(x, qmin, otherParams, tau, rho){
    n.r <- length(rho)
    n.x <- length(x)
    
    t1 <- tau %x% t(x)
    tmp <- (qmin - t1)/(1-rho)
    minval <- apply(tmp,2,min)

    degf <- otherParams["degf"]
    mu <- otherParams["mu"]
    varia <- otherParams["varia"]

    temp.q<-(minval - mu)/sqrt(varia)*sqrt(2*degf) + degf

    re<-pchisq(temp.q ,df=degf) * dchisq(x,df=1)

    return(re)
}
