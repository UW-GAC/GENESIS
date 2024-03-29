
## function that gets an n \times p matrix of p genotypes of n individuals, and a null model, and tests the genotypes associations with the outcomes.
## Genetic data are always assumed complete.
## Types of tests:
## Variant set: SKAT, burden, SKAT-O. Multiple types of p-values. Default: Davies with Kuonen if does not converge.


testVariantSet <- function( nullmod, G, weights,
                            test = c("Burden", "SKAT", "fastSKAT", "SMMAT", "fastSMMAT", "SKATO", "BinomiRare", "CMP"),
                            # burden.test = c("Score"),
                            neig = 200, ntrace = 500,
                            rho = seq(from = 0, to = 1, by = 0.1)){
                           # pval.method = c("davies", "kuonen", "liu"),
                           # return.scores = FALSE, return.scores.cov = FALSE){

    test <- match.arg(test)
    # burden.test <- match.arg(burden.test)
    # pval.method <- match.arg(pval.method)

    G <- .genoAsMatrix(nullmod, G)

    if (test == "Burden") {
        out <- .testVariantSetBurden(nullmod, G, weights, burden.test = "Score")
    }
    if (test == "BinomiRare") {
        out <- .testVariantSetBurden(nullmod, G, weights, burden.test = "BinomiRare")
    }
    if (test == "CMP") {
        out <- .testVariantSetBurden(nullmod, G, weights, burden.test = "CMP")
    }
    if (test == "SKAT") {
        out <- .testVariantSetSKAT(nullmod, G, weights, neig = Inf, ntrace = Inf)
                                   # return.scores, return.scores.cov)
    }
    if(test == "fastSKAT"){
        out <- .testVariantSetSKAT(nullmod, G, weights, neig, ntrace)
    }
    if (test == "SMMAT") {
        out <- .testVariantSetSMMAT(nullmod, G, weights, neig = Inf, ntrace = Inf)
    }
    if(test == "fastSMMAT"){
        out <- .testVariantSetSMMAT(nullmod, G, weights, neig, ntrace)
    }
    if(test == "SKATO"){
        out <- .testVariantSetSKATO(nullmod, G, weights, rho)
    }
    return(out)
}



## create the burden score, than calls the appropriate single variant test function.
## can easily implement GxE interaction with the burden score... later!
.testVariantSetBurden <- function(nullmod, G, weights, burden.test){
    # multiply G by weights and compute burden
    if(is(G, "Matrix")){
        burden <- rowSums(G %*% Diagonal(x = weights))
    }else{
        burden <- colSums(t(G) * weights)
    }

    # adjust burden for covariates and random effects
    if (burden.test == "Score") {
        Gtilde <- calcGtilde(nullmod, burden)
        out <- .testGenoSingleVarScore(Gtilde, G = burden, resid = nullmod$fit$resid.PY, RSS0 = nullmod$RSS0)
    }
    # if (burden.test == "Wald"){
    #     out <- .testGenoSingleVarWald(Gtilde, Ytilde = nullmod$Ytilde,
    #                                   n = length(nullmod$Ytilde), k = ncol(nullmod$model.matrix))
    # }

    if (burden.test == "BinomiRare"){
        ## if this is a mixed model, used conditional probabilities
        if (nullmod$model$family$mixedmodel) {
            phat <- expit(nullmod$fit$linear.predictor)
        } else {
            phat <- nullmod$fit$fitted.values
        }
        out <- .testGenoSingleVarBR(nullmod$fit$outcome, probs=phat, G=matrix(burden))
    }

    if (burden.test == "CMP"){
        if (nullmod$model$family$mixedmodel) {

            phat <- expit(nullmod$fit$linear.predictor)
            out <- .testGenoSingleVarCMP(nullmod$fit$outcome, probs=phat, G=matrix(burden))
        } else{ ## not a mixed model
            phat <- nullmod$fit$fitted.values
            out <- .testGenoSingleVarBR(nullmod$fit$outcome, probs=phat, G=matrix(burden))
        }
    }
    return(out)
}


## new function that runs both SKAT and fastSKAT
.testVariantSetSKAT <- function(nullmod, G, weights, neig = 200, ntrace = 500, verbose = FALSE){
    # multiply G by weights
    if(is(G, "Matrix")){
        G <- G %*% Diagonal(x = weights)
    }else{
        G <- t(t(G) * weights)
    }

    # scores
    U <- as.vector(crossprod(G, nullmod$fit$resid.PY)) # WGPY
    # SKAT test statistic
    Q <- sum(U^2)

    # adjust G for covariates and random effects
    G <- calcGtilde(nullmod, G) # P^{1/2}GW

    # compute the p-value
    out <- .calcPvalVCTest(Q = Q, G = G, neig = neig, ntrace = ntrace, verbose = verbose)

    return(list(Q = Q, pval = out$pval, err = out$err, pval.method = out$pval.method))
}

## function for SMMAT and fastSMMAT
.testVariantSetSMMAT <- function(nullmod, G, weights, neig = 200, ntrace = 500, verbose = FALSE) {
    # multiply G by weights
    if(is(G, "Matrix")){
        G <- G %*% Diagonal(x = weights)
    }else{
        G <- t(t(G) * weights)
    }

    # scores
    U <- as.vector(crossprod(G, nullmod$fit$resid.PY)) # WGPY
    U.sum <- sum(U) # 1WGPY

    # adjust G for covariates and random effects
    G <- calcGtilde(nullmod, G) # P^{1/2}GW

    # compute burden p-value
    G.rowSums <- rowSums(G) # P^{1/2}GW1
    GG1 <- crossprod(G, G.rowSums) # WGPGW1  # O(mn)
    V.sum <- sum(GG1) # 1WGPGW1
    burden.stat = U.sum / sqrt(V.sum)
    burden.pval <- .pchisq_filter_extreme(burden.stat^2, df=1, lower.tail=FALSE)

    # adjust U and G for burden
    U <- U - GG1*U.sum/V.sum # WGPY - WGPGW1 * 1WGPY/(1WGPGW1)
    G <- G - tcrossprod(G.rowSums, GG1)/V.sum # O(mn)

    # SMMAT test statistic
    Q <- sum(U^2)

    ### alternative to part above; seems to be slower from testing; this is how presented in SMMAT paper ###
    # V <- crossprod(G) # WGPGW  # O(m^2n)
    # GG1 <- rowSums(V) # WGPGW1
    # # denominator for burden
    # V.sum <- sum(GG1) # 1WGPGW1
    # # burden p-value
    # burden.pval <- .pchisq_filter_extreme(U.sum^2/V.sum, df=1, lower.tail=FALSE)
    # # adjust for burden
    # U <- U - GG1*U.sum/V.sum
    # V <- V - tcrossprod(GG1)/V.sum  # O(m^2)

    # compute the p-value for the "adjusted SKAT" part
    out <- .calcPvalVCTest(Q = Q, G = G, neig = neig, ntrace = ntrace, verbose = verbose)
    theta.pval <- out$pval
    err <- out$err

    # Fisher's method to combine p-values
    smmat.pval <- tryCatch(.pchisq_filter_extreme(-2*log(burden.pval)-2*log(theta.pval), df=4, lower.tail = FALSE), error = function(e) { NA })
    if(is.na(smmat.pval)) {
        err <- 1
        smmat.pval <- burden.pval
    }
    return(list(
      # burden test output.
      Score_burden = U.sum,
      Score.SE_burden = sqrt(V.sum),
      Stat_burden = burden.stat,
      pval_burden = burden.pval,
      # SKAT output.
      # SMMAT output.
      Q_theta = Q,
      pval_theta = theta.pval,
      pval_SMMAT = smmat.pval,
      err = err,
      pval_theta.method = out$pval.method
    ))
}


.regular <- function(Q, V, ncolG) {
    if(ncolG == 1){
        pv <- list(pval = .pchisq_filter_extreme(as.numeric(Q/V), df=1, lower.tail=FALSE), method = "integration")

    }else{
        lambda <- eigen(V, only.values = TRUE, symmetric=TRUE)$values
        # lambda <- lambda[lambda > 0]
        pv <- .pchisqsum(x = Q, df = rep(1, length(lambda)), a = lambda)
        # pv <- tryCatch({
        #             list(pval = .pchisqsum(x = Q, df = rep(1, length(lambda)), a = lambda, method = "integration"),
        #                  method = "integration")
        #         }, warning = function(w){
        #             list(pval = pchisqsum(x = Q, df = rep(1, length(lambda)), a = lambda, method = "saddlepoint"),
        #                  method = "saddlepoint")
        #         }, error = function(e){
        #             list(pval = NA_real_,
        #                  method = NA_character_)
        #         })
    }
    pv$err <- ifelse(is.na(pv$pval), 1, 0)
    return(pv)
}

.fastH <- function(Q, V, neig) {
    pv <- list(pval = NA_real_, method = NA_character_)
    pval.try = 0
    while(is.na(pv$pval) & pval.try < 10){
        pv <- tryCatch( pchisqsum_ssvd(x = Q, M = as.matrix(V), n = neig, p = 10, q = 1),
                       error = function(e){ list(pval = NA_real_, method = "error") } )
        pval.try <- pval.try + 1
    }
    pv$method <- paste0('ssvd_', pv$method)
    if(is.na(pv$pval)){
        err <- 1
    }else if(pval.try > 1){
        err <- 2
    }else{
        err <- 0
    }
    pv[["err"]] <- err
    return(pv)
}

.fastG <- function(Q, G, neig, ntrace) {
    pv <- list(pval = NA_real_, method = NA_character_)
    pval.try = 0
    while(is.na(pv$pval) & pval.try < 10){
        pv <- tryCatch( pchisqsum_rsvd(x = Q, M = as.matrix(G), n = neig, p = 10, q = 3, tr2.sample.size = ntrace),
                       error = function(e){ list(pval = NA_real_, method = "error") }   )
        pval.try <- pval.try + 1
    }
    pv$method <- paste0('rsvd_', pv$method)
    if(is.na(pv$pval)){
        err <- 1
    }else if(pval.try > 1){
        err <- 2
    }else{
        err <- 0
    }
    pv[["err"]] <- err
    return(pv)
}

.calcPvalVCTest <- function(Q, G, neig, ntrace, verbose){
    if(!requireNamespace("survey")) stop("package 'survey' must be installed to calculate p-values for SKAT or SMMAT")
    if(!requireNamespace("CompQuadForm")) stop("package 'CompQuadForm' must be installed to calculate p-values for SKAT or SMMAT")

    ncolG <- ncol(G) # number of snps
    nrowG <- nrow(G) # number of samples
    if (verbose) message('nsamp = ', nrowG, '; nsnp = ', ncolG)

    if(min(ncolG, nrowG) < 6000 + 20*neig){
        if(ncolG <= nrowG){
            V <- crossprod(G) # WGPGW
        }else{
            V <- tcrossprod(G) # same eigenspace but smaller matrix
        }
        if(mean(abs(V)) < sqrt(.Machine$double.eps)){
            return(list(pval = NA_real_, pval.method = NA_character_, err = 1))
        }

        if(min(ncolG, nrowG) < 2*neig){
            # use "regular" method
            pv <- .regular(Q, V, ncolG)

        }else{
            # use "fast H" method
            if (verbose) message("using method fast_H")
            pv <- .fastH(Q, V, neig)
        }

    }else{
        # use "fast G" method
        if (verbose) message("using method fast_G")
        pv <- .fastG(Q, G, neig, ntrace)
    }

    return(list(pval = pv$pval, pval.method = pv$method, err = pv$err))
}

.pchisqsum <- function(x, df, a){

    ## check for bad.df
    ## can happen with randomised trace estimator if most remaining singular values are very small
    ## leads to unreliable p-values
    if(any(df < 1)){
        stop("Negative/fractional df")
    }

    df<-round(df)

    ## try integration
    f <- suppressWarnings(CompQuadForm::davies(x, a, df, acc = 1e-9))
    if((f$ifault > 0) | (f$Qq < 1e3*.Machine$double.eps) | (f$Qq > 1)){
        ## try saddlepoint
        pval <- survey:::saddle(x, rep(a, df))
        method <- "saddlepoint"
    }else{
        pval <- f$Qq
        method <- "integration"
    }

    return(list(pval = pval, method = method))
}

# .calcPval <- function(Q, lambda, pval.method) {
#     if(!requireNamespace("survey")) stop("package 'survey' must be installed to calculate p-values for SKAT")
#     if(!requireNamespace("CompQuadForm")) stop("package 'CompQuadForm' must be installed to calculate p-values for SKAT")

#     err <- 0
#     if(pval.method == "kuonen"){
#         pval <- survey:::saddle(x = Q, lambda = lambda)
#         err <- ifelse(is.na(pval), 1, 0)

#     }else if(pval.method == "davies"){
#         tmp <- suppressWarnings(CompQuadForm::davies(q = Q, lambda = lambda, acc = 1e-06))
#         pval <- tmp$Qq
#         if((tmp$ifault > 0) | (pval <= 0) | (pval >= 1)) {
#             pval <- survey:::saddle(x = Q, lambda = lambda)
#         }
#         err <- ifelse(is.na(pval), 1, 0)

#     }else if(pval.method == "liu"){
#         pval <- CompQuadForm::liu(q = Q, lambda = lambda)
#         err <- 0
#     }

#     if(err > 0){
#         pval <- CompQuadForm::liu(q = Q, lambda = lambda)
#     }

#     return(c(pval=pval, err=err))
# }


## new function just for SKAT-O
.testVariantSetSKATO <- function(nullmod, G, weights, rho = 0){
#                                 # return.scores = FALSE, return.scores.cov = FALSE){
    # scores
    scores <- as.vector(crossprod(G, nullmod$fit$resid.PY))

    # adjust G for covariates and random effects
    geno.adj <- calcGtilde(nullmod, G)

    # covariance of scores
    V <- crossprod(geno.adj)

    # vectors to hold output
    nrho <- length(rho)
    out.Q <- rep(NA, nrho); names(out.Q) <- paste("Q", rho, sep="_")
    out.pval <- rep(NA, nrho); names(out.pval) <- paste("pval", rho, sep="_")
    out.err <- rep(NA, nrho); names(out.err) <- paste("err", rho, sep="_")
    lambdas <- vector("list", nrho)

    # get p-value for each choice of rho
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

        # p value calculation
        if(length(scores) == 1){
            lambdas[[i]] <- as.numeric(distMat)
            pval <- .pchisq_filter_extreme(as.numeric(Q/distMat), df=1, lower.tail=FALSE)
            err <- ifelse(is.na(pval), 1, 0)
        }else{
            lambda <- eigen(distMat, only.values = TRUE, symmetric=TRUE)$values
            # lambda <- lambda[lambda > 0]
            lambdas[[i]] <- lambda
            pv <- .pchisqsum(x = Q, df = rep(1, length(lambda)), a = lambda)
            pval <- pv$pval
            err <- ifelse(is.na(pval), 1, 0)
        }

        # update results
        out.Q[i] <- Q
        out.pval[i] <- pval
        out.err[i] <- err
    }
    out <- as.list(c(out.Q, out.pval, out.err))

    # get SKAT-O p-value
    if(length(scores) == 1){
        # pvalue is the same for all rhos
        out2 <- list(min.pval=out.pval[1], opt.rho=NA, pval_SKATO=out.pval[1])

    }else{
        # find the minimum p-value
        minp <- min(out.pval, na.rm=TRUE)
        opt.rho <- if (!is.na(minp)) rho[which.min(out.pval)] else NA
        out2 <- list(min.pval=minp, opt.rho=opt.rho)

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

    # return results
    return(out)
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

    re<-.pchisq_filter_extreme(temp.q ,df=degf) * dchisq(x,df=1)

    return(re)
}
