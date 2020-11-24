
## function that gets an n\times p matrix of p genotypes of n individuals, and a null model, and tests the genotypes associations with the outcomes.
## Genetic data are always assumed complete.
## Types of tests:
## Single variant: Score, Score.SPA, BinomiRare, interaction.
## Variant set: SKAT, burden, SKAT-O. Multiple types of p-values. Default: Davis with Koenen if does not converge.


# E an environmental variable for optional GxE interaction analysis.
testGenoSingleVar <- function(nullmod, G, E = NULL, test = c("Score", "Score.SPA", "BinomiRare", "CMP"),
                              recalc.pval.thresh = 1, GxE.return.cov = FALSE){
    test <- match.arg(test)
    calc.score <- test %in% c("Score", "Score.SPA") | (recalc.pval.thresh < 1)

    if (isNullModelSmall(nullmod) && (calc.score || !is.null(E))) {
        stop("small null model cannot be used with options provided")
    }

    G <- .genoAsMatrix(nullmod, G)

    # checks on test
    if (!is.null(E)){
        #message("Performing GxE test")
        res <- .testGenoSingleVarWaldGxE(nullmod, G, E, GxE.return.cov.mat=GxE.return.cov)
        return(res)
    }

    if(test == "Score.SPA" & nullmod$model$family$family != "binomial"){
        test <- "Score"
        message("Saddlepoint approximation (SPA) can only be used for binomial family; using Score test instead.")
    }

    # run the test
    if(calc.score){
        Gtilde <- calcGtilde(nullmod, G)
        res <- .testGenoSingleVarScore(Gtilde, G, nullmod$fit$resid.PY, nullmod$RSS0)
    }

    if(test == "Score.SPA"){
        # saddle point approximation
        res <- SPA_pval(score.result = res, nullmod = nullmod, G = G, pval.thresh = recalc.pval.thresh)
    }

    if (test == "BinomiRare"){
      if (nullmod$model$family$family != "binomial") stop("BinomiRare can only be used for binomial family.")

      if (nullmod$model$family$mixedmodel) { ## if this is a mixed model, use conditional probabilities $
        phat <- expit(nullmod$fit$workingY - nullmod$fit$resid.conditional)
      } else{ ## not a mixed model
        phat <- nullmod$fit$fitted.values
      }
      score.pval <- if(calc.score) res$Score.pval else NULL
      res <- .testGenoSingleVarBR(nullmod$fit$outcome, probs=phat, G, score.pval=score.pval, pval.thresh=recalc.pval.thresh)
    }

    if (test == "CMP"){
      score.pval <- if(calc.score) res$Score.pval else NULL
      if (nullmod$model$family$mixedmodel) { ## if this is a mixed model, use conditional probabilities.
        phat <- expit(nullmod$fit$workingY - nullmod$fit$resid.conditional)
        res <- .testGenoSingleVarCMP(nullmod$fit$outcome, probs=phat, G, score.pval=score.pval, pval.thresh=recalc.pval.thresh)
      } else{ ## not a mixed model
        phat <- nullmod$fit$fitted.values
        res <- .testGenoSingleVarBR(nullmod$fit$outcome, probs=phat, G,  score.pval=score.pval, pval.thresh=recalc.pval.thresh)
      }
    }

    return(res)
}



## this function currently assumes that the alt allele is the minor allele. So either G
## needs to be such that alt allele is minor allele, or the function checks for it, or a vector of
## indicators or of frequencies would be provided.
.testGenoSingleVarBR <- function(D, probs, G, score.pval=NULL, pval.thresh = 0.05){
    if (!requireNamespace("poibin")) stop("package 'poibin' must be installed for the BinomiRare test")
    cols <- c("n.carrier", "n.D.carrier", "expected.n.D.carrier", "pval", "mid.pval")
    res <- matrix(NA, nrow = ncol(G), ncol = length(cols), dimnames = list(NULL, cols))

    for (i in seq(ncol(G))){
        if (sd(G[,i])==0){
            next
        }
        carrier.inds <- which(G[,i] > 0)
        res[i, "n.carrier"] <- length(carrier.inds)
        cur.prob.vec <- probs[carrier.inds]
        res[i, "expected.n.D.carrier"] <- sum(cur.prob.vec)
        res[i, "n.D.carrier"] <- sum(D[carrier.inds])
        if (!is.null(score.pval) && score.pval[i] > pval.thresh){
            res[i, c("pval", "mid.pval")] <- score.pval[i]
        } else{
            res[i, c("pval", "mid.pval")] <- .poibinP(n.carrier = length(carrier.inds), n.D.carrier = sum(D[carrier.inds]), prob.vec = cur.prob.vec)
        }
    }
    res <- as.data.frame(res)
    return(res)
}


.poibinP <- function(n.carrier, n.D.carrier, prob.vec){
    stopifnot(n.D.carrier <= n.carrier, length(prob.vec) == n.carrier)
    d.poibin <- poibin::dpoibin(0:n.carrier, prob.vec)
    prob.cur <- d.poibin[n.D.carrier + 1]
    pval <- prob.cur + sum(d.poibin[d.poibin < prob.cur])
    mid.pval <- 0.5*prob.cur + sum(d.poibin[d.poibin < prob.cur])

    pvals <- c(pval=pval, mid.pval=mid.pval)
    return(pvals)
}



.testGenoSingleVarScore <- function(Gtilde, G, resid, RSS0){
    GPG <- colSums(Gtilde^2) # vector of G^T P G (for each SNP)
    score.SE <- sqrt(GPG)
    score <- as.vector(crossprod(G, resid)) # G^T P Y
    Stat <- score/score.SE

    res <- data.frame(Score = score, Score.SE = score.SE, Score.Stat = Stat,
                      Score.pval = pchisq(Stat^2, df = 1, lower.tail = FALSE),
                      Est = score/GPG, Est.SE = 1/score.SE,
                      PVE = (Stat^2)/RSS0) # RSS0 = (n-k) when gaussian; not when binary

    return(res)
}



# .testGenoSingleVarWald <- function(Gtilde, Ytilde, n, k){
#     GPG <- colSums(Gtilde^2) # vector of G^T P G (for each SNP)
#     GPY <- as.vector(crossprod(Gtilde, Ytilde)) # vector of G^T P Y (for each SNP)
#     beta <- GPY/GPG
#     sY2 <- sum(Ytilde^2)
#     RSS <- as.numeric((sY2 - GPY * beta)/(n - k - 1))
#     Vbeta <- RSS/GPG
#     Stat <- beta/sqrt(Vbeta)
#     res <- data.frame(Est = beta, Est.SE = sqrt(Vbeta), Wald.Stat = Stat,
#                       Wald.pval = pchisq(Stat^2, df = 1, lower.tail = FALSE))
#     return(res)
# }


.testGenoSingleVarWaldGxE <- function(nullmod, G, E, GxE.return.cov.mat = FALSE){

    E <- as.matrix(E)
    p <- ncol(G)
    v <- ncol(E) + 1
    n <- nrow(nullmod$fit)
    k <- ncol(nullmod$model.matrix)
    sY2 <- as.numeric(crossprod(nullmod$fit$resid.cholesky))

    if (GxE.return.cov.mat) {
        res.Vbetas <- vector(mode = "list", length = p)
    }

    intE <- cbind(1, E) # add intercept the "Environmental" variable E.
    if (is(G, "Matrix")) intE <- Matrix(intE)

    var.names <- c("G", paste("G", colnames(E), sep = ":"))

    res <- matrix(NA, nrow = p, ncol = length(var.names)*2 + 2,
                  dimnames = list(NULL,
                                  c(paste0("Est.", var.names), paste0("SE.", var.names), "GxE.Stat", "Joint.Stat" ) ))

    for (g in 1:p) {
        Gtilde <- calcGtilde(nullmod, G[, g] * intE)
        GPG <- crossprod(Gtilde)
        GPGinv <- tryCatch(chol2inv(chol(GPG)), error = function(e) {TRUE}) # this is inverse A matrix of sandwich
        # check that the error function above hasn't been called (which returns TRUE instead of the inverse matrix)
        if (is.logical(GPGinv)) next

        GPY <- crossprod(Gtilde, nullmod$fit$resid.cholesky)
        betas <- crossprod(GPGinv, GPY)
        res[g, grep("^Est\\.G", colnames(res))] <- as.vector(betas)

        RSS <- as.numeric((sY2 - crossprod(GPY, betas))/(n - k - v))
        Vbetas <- GPGinv * RSS

        if (GxE.return.cov.mat) {
            res.Vbetas[[g]] <- Vbetas
        }

        res[g, grep("^SE\\.G", colnames(res))] <- sqrt(diag(Vbetas))

        res[g, "GxE.Stat"] <- tryCatch(sqrt(as.vector(crossprod(betas[-1],
                                                 crossprod(chol2inv(chol(Vbetas[-1, -1])),
                                                           betas[-1])))),
                                       error = function(e) { NA })

        res[g, "Joint.Stat"] <- tryCatch(sqrt(as.vector(crossprod(betas,
                                                   crossprod(GPG, betas))/RSS)),
                                         error = function(e) { NA })
    }

    res <- as.data.frame(res)
    res$GxE.pval <- pchisq((res$GxE.Stat)^2, df = (v - 1), lower.tail = FALSE)
    res$Joint.pval <- pchisq((res$Joint.Stat)^2, df = v, lower.tail = FALSE)

    if (GxE.return.cov.mat) {
        return(list(res = res, GxEcovMatList = res.Vbetas))
    } else {
        return(res)
    }
}



.testGenoSingleVarCMP <- function(D, probs, G, score.pval=NULL, pval.thresh = 0.05){
    if (!requireNamespace("COMPoissonReg")) stop("package 'COMPoissonReg' must be installed for the CMP test")
    cols <- c("ncoln.carrier", "n.D.carrier", "expected.n.D.carrier", "pval", "mid.pval")
    res <- matrix(NA, nrow = ncol(G), ncol = length(cols), dimnames = list(NULL, cols))

    for (i in seq(ncol(G))){
        if (sd(G[,i])==0){
            next
        }
        carrier.inds <- which(G[,i] > 0)
        phat <- probs[carrier.inds]
        ncar <- length(phat)
        sum.d <- sum(D[carrier.inds])

        res[i, "n.carrier"] <- length(carrier.inds)
        res[i, "n.D.carrier"] <- sum.d
        cur.prob.vec <- probs[carrier.inds]
        res[i, "expected.n.D.carrier"] <- sum(cur.prob.vec)
        if (!is.null(score.pval) && score.pval[i] > pval.thresh){
            res[i, c("pval", "mid.pval")] <- score.pval[i]
        } else {
            if (ncar == 1) {
                res[i, c("pval", "mid.pval")] <- ifelse(sum.d == 1, phat, 1-phat)
                next
            }

            mu1.analytic <- sum(phat)
            var.analytic <- sum(phat*(1-phat))
            nuhat <- mu1.analytic/var.analytic
            lamhat <- mu1.analytic^nuhat
            res[i, c("pval", "mid.pval")] <- .calc_cmp_pval(ncar, sum.d, lamhat, nuhat)
        }

    }
    res <- as.data.frame(res)
    return(res)
}


## compute p-value for CMP test based on estimated lambda, nu, number of carriers, and nubmber of diseased carriers.
.calc_cmp_pval <- function(ncar, sum.d, lamhat, nuhat){ #, midp.type = "both"){
    prob.cur <- COMPoissonReg::dcmp(sum.d + 1, lamhat, nuhat)
    d.cmp <- COMPoissonReg::dcmp(0:ncar, lamhat, nuhat)
    pval <- prob.cur + sum(d.cmp[d.cmp < prob.cur])
    mid.pval <- pval-prob.cur/2
    pvals <- c(pval=pval, mid.pval=mid.pval)
    return(pval)
}



## G is an n by v matrix of 2 or more columns, all representing alleles of the same (multi-allelic) variant.
.testSingleVarMultAlleles <- function(Gtilde, Ytilde, n, k){
    v <- ncol(Gtilde)

    var.names <- colnames(Gtilde)

    res <- matrix(NA, nrow = 1, ncol = length(var.names)*2 + 2,
                  dimnames = list(NULL,
                                  c(paste0("Est.", var.names), paste0("SE.", var.names), "Joint.Stat", "Joint.Pval" ) ))


    GPG <- crossprod(Gtilde)
    GPGinv <- tryCatch(chol2inv(chol(GPG)), error = function(e) {TRUE})

    if (is.logical(GPGinv)) return(list(res = res, allelesCovMat = NA))

    GPY <- crossprod(Gtilde, Ytilde)
    betas <- crossprod(GPGinv, GPY) ## effect estimates of the various alleles
    res[1, grep("^Est\\.G", colnames(res))] <- betas

    sY2 <- sum(Ytilde^2)
    RSS <- as.numeric((sY2 - crossprod(GPY, betas))/(n - k - v))
    Vbetas <- GPGinv * RSS

    res[1, grep("^SE\\.G", colnames(res))] <- sqrt(diag(Vbetas))

    res[1, "Joint.Stat"] <- tryCatch(crossprod(betas,
                                               crossprod(GPG, betas))/RSS,
                                     error = function(e) { NA })

    res[,"Joint.pval"] <- pchisq(res[,"Joint.Stat"], df = v, lower.tail = FALSE)

    res <- as.data.frame(res)
    return(list(res = res, allelesCovMat = Vbetas))
}


## include expit function (not included in base R)
expit <- function(x){
    ex <- exp(x)
    ex/(1+ex)
}
