context("check variant set association tests")

.testNullmod <- function(n, MM=FALSE, binary=FALSE) {
    X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
    if (!binary) {
	y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))
        family <- "gaussian"
    }  else {
        y <- rbinom(n, size = 1, prob = 0.4)
        family <- "binomial"
    }
    
    group.idx <- list(G1 = c(1:(n/2)), G2 = c((n/2 + 1):n))

    if (MM) {
        cor.mat <- matrix(rnorm(n*n, sd = 0.05),n,n)
        cor.mat <- crossprod(cor.mat)
        covMatList <- list(A = cor.mat)
    } else {
        covMatList <- NULL
    }
    
    fitNullMod(y, X, covMatList = covMatList, group.idx = group.idx, family=family, verbose=FALSE)
}


test_that("SKAT with rho=1 matches burden", {
	n <- 100

	### create a matrix of genetic variants to test.
	geno <- matrix(rbinom(200*n, size = 2, prob = 0.2), nrow = n, ncol = 200)
        weights <- rep(1, ncol(geno))

        ## mixed model
	nullmod <- .testNullmod(n, MM=TRUE)
        skat <- .testVariantSetSKAT(nullmod, geno, weights, rho=1, pval.method="davies")
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_true(abs(skat$pval_1 - burden$Score.pval) < 0.01)
        
        ## basic
	nullmod <- .testNullmod(n, MM=FALSE)
        skat <- .testVariantSetSKAT(nullmod, geno, weights, rho=1, pval.method="davies")
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_true(abs(skat$pval_1 - burden$Score.pval) < 0.01)
        
        ## mixed model - binary
	nullmod <- .testNullmod(n, MM=TRUE, binary=TRUE)
        skat <- .testVariantSetSKAT(nullmod, geno, weights, rho=1, pval.method="davies")
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_true(abs(skat$pval_1 - burden$Score.pval) < 0.01)
        
        ## basic - binary
	nullmod <- .testNullmod(n, MM=FALSE, binary=TRUE)
        skat <- .testVariantSetSKAT(nullmod, geno, weights, rho=1, pval.method="davies")
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_true(abs(skat$pval_1 - burden$Score.pval) < 0.01)
})


test_that("hybrid test matches burden and skat", {
	n <- 100

	### create a matrix of genetic variants to test.
	geno <- matrix(rbinom(200*n, size = 2, prob = 0.2), nrow = n, ncol = 200)
        weights <- rep(1, ncol(geno))

        ## mixed model
	nullmod <- .testNullmod(n, MM=TRUE)
        hybrid <- .testVariantSetSMMAT(nullmod, geno, weights, pval.method="davies")
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(hybrid$pval_burden, burden$Score.pval)
        
        ## basic
	nullmod <- .testNullmod(n, MM=FALSE)
        hybrid <- .testVariantSetSMMAT(nullmod, geno, weights, pval.method="davies")
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(hybrid$pval_burden, burden$Score.pval)
        
        ## mixed model - binary
	nullmod <- .testNullmod(n, MM=TRUE, binary=TRUE)
        hybrid <- .testVariantSetSMMAT(nullmod, geno, weights, pval.method="davies")
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(hybrid$pval_burden, burden$Score.pval)
        
        ## basic - binary
	nullmod <- .testNullmod(n, MM=FALSE, binary=TRUE)
        hybrid <- .testVariantSetSMMAT(nullmod, geno, weights, pval.method="davies")
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(hybrid$pval_burden, burden$Score.pval)
})
