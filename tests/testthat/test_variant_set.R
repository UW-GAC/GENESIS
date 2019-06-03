context("check variant set association tests")

test_that("SKAT-O with rho=1 matches burden", {
	n <- 100

	### create a matrix of genetic variants to test.
	geno <- .testGenoMatrix(n=n)
        weights <- rep(1, ncol(geno))

        ## mixed model
	nullmod <- .testNullmod(n, MM=TRUE)
        skato <- .testVariantSetSKATO(nullmod, geno, weights, rho=1)
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(skato$pval_1, burden$Score.pval, tolerance=0.01)
        
        ## basic
	nullmod <- .testNullmod(n, MM=FALSE)
        skato <- .testVariantSetSKATO(nullmod, geno, weights, rho=1)
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(skato$pval_1, burden$Score.pval, tolerance=0.01)
        
        ## mixed model - binary
	nullmod <- .testNullmod(n, MM=TRUE, binary=TRUE)
        skato <- .testVariantSetSKATO(nullmod, geno, weights, rho=1)
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(skato$pval_1, burden$Score.pval, tolerance=0.01)
        
        ## basic - binary
	nullmod <- .testNullmod(n, MM=FALSE, binary=TRUE)
        skato <- .testVariantSetSKATO(nullmod, geno, weights, rho=1)
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(skato$pval_1, burden$Score.pval, tolerance=0.01)
})


test_that("SMMAT matches burden", {
	n <- 100

	### create a matrix of genetic variants to test.
	geno <- .testGenoMatrix(n=n)
        weights <- rep(1, ncol(geno))

        ## mixed model
	nullmod <- .testNullmod(n, MM=TRUE)
        smmat <- .testVariantSetSMMAT(nullmod, geno, weights)
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(smmat$pval_burden, burden$Score.pval)
        
        ## basic
	nullmod <- .testNullmod(n, MM=FALSE)
        smmat <- .testVariantSetSMMAT(nullmod, geno, weights)
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(smmat$pval_burden, burden$Score.pval)
        
        ## mixed model - binary
	nullmod <- .testNullmod(n, MM=TRUE, binary=TRUE)
        smmat <- .testVariantSetSMMAT(nullmod, geno, weights)
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(smmat$pval_burden, burden$Score.pval)
        
        ## basic - binary
	nullmod <- .testNullmod(n, MM=FALSE, binary=TRUE)
        smmat <- .testVariantSetSMMAT(nullmod, geno, weights)
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(smmat$pval_burden, burden$Score.pval)
})


test_that("SKATO single rho", {
	n <- 100

	### create a matrix of genetic variants to test.
	geno <- .testGenoMatrix(n=n)
        weights <- rep(1, ncol(geno))

        ## mixed model
	nullmod <- .testNullmod(n, MM=TRUE)
        skato <- .testVariantSetSKATO(nullmod, geno, weights, rho=0)
        expect_equal(skato$pval_0, skato$min.pval)
        expect_equal(skato$pval_0, skato$pval_SKATO, tolerance=0.001)
})


test_that("SKATO multiple rho", {
	n <- 100

	### create a matrix of genetic variants to test.
	geno <- .testGenoMatrix(n=n)
        weights <- rep(1, ncol(geno))

        ## mixed model
	nullmod <- .testNullmod(n, MM=TRUE)
        skato <- .testVariantSetSKATO(nullmod, geno, weights, rho=c(0,0.5,1))
        expect_true(all(paste0("pval_", c(0,0.5,1)) %in% names(skato)))
})


test_that("SKAT matches SKAT-O with rho=0", {
	n <- 100

	### create a matrix of genetic variants to test.
	geno <- .testGenoMatrix(n=n)
        weights <- rep(1, ncol(geno))

        ## mixed model
	nullmod <- .testNullmod(n, MM=TRUE)
        skat <- .testVariantSetSKAT(nullmod, geno, weights)
        skato <- .testVariantSetSKATO(nullmod, geno, weights, rho=0)
        expect_equal(skat$pval, skato$pval_0)
        
        ## basic
	nullmod <- .testNullmod(n, MM=FALSE)
        skat <- .testVariantSetSKAT(nullmod, geno, weights)
        skato <- .testVariantSetSKATO(nullmod, geno, weights, rho=0)
        expect_equal(skat$pval, skato$pval_0)
        
        ## mixed model - binary
	nullmod <- .testNullmod(n, MM=TRUE, binary=TRUE)
        skat <- .testVariantSetSKAT(nullmod, geno, weights)
        skato <- .testVariantSetSKATO(nullmod, geno, weights, rho=0)
        expect_equal(skat$pval, skato$pval_0)
        
        ## basic - binary
	nullmod <- .testNullmod(n, MM=FALSE, binary=TRUE)
        skat <- .testVariantSetSKAT(nullmod, geno, weights)
        skato <- .testVariantSetSKATO(nullmod, geno, weights, rho=0)
        expect_equal(skat$pval, skato$pval_0)
})


test_that("fastSKAT(H) matches SKAT-O with rho=0", {
	n <- 250

	### create a matrix of genetic variants to test.
	geno <- .testGenoMatrix(n=n, nsnp=300)
        weights <- rep(1, ncol(geno))

        ## mixed model
	nullmod <- .testNullmod(n, MM=TRUE)
        expect_message(skat <- .testVariantSetSKAT(nullmod, geno, weights, neig=100, verbose=TRUE), "fast_H")
        skato <- .testVariantSetSKATO(nullmod, geno, weights, rho=0)
        expect_equal(skat$pval, skato$pval_0, tolerance=0.01)
        
        ## basic
	nullmod <- .testNullmod(n, MM=FALSE)
        expect_message(skat <- .testVariantSetSKAT(nullmod, geno, weights, neig=100, verbose=TRUE), "fast_H")
        skato <- .testVariantSetSKATO(nullmod, geno, weights, rho=0)
        expect_equal(skat$pval, skato$pval_0, tolerance=0.01)
        
        ## mixed model - binary
	nullmod <- .testNullmod(n, MM=TRUE, binary=TRUE)
        expect_message(skat <- .testVariantSetSKAT(nullmod, geno, weights, neig=100, verbose=TRUE), "fast_H")
        skato <- .testVariantSetSKATO(nullmod, geno, weights, rho=0)
        expect_equal(skat$pval, skato$pval_0, tolerance=0.01)
        
        ## basic - binary
	nullmod <- .testNullmod(n, MM=FALSE, binary=TRUE)
        expect_message(skat <- .testVariantSetSKAT(nullmod, geno, weights, neig=100, verbose=TRUE), "fast_H")
        skato <- .testVariantSetSKATO(nullmod, geno, weights, rho=0)
        expect_equal(skat$pval, skato$pval_0, tolerance=0.01)
})


test_that("fastSKAT matches SKAT(regular)", {
	n <- 500
	nullmod <- .testNullmod(n, MM=TRUE)
	G <- .testGenoMatrix(n=n, nsnp=300)
        U <- as.vector(crossprod(G, nullmod$resid))
        Q <- sum(U^2)
        G <- calcGtilde(nullmod, G)
        V <- tcrossprod(G)

        ## mixed model
        reg <- .regular(Q, V, ncol(G))
        fastH <- .fastH(Q, V, neig=200)
        expect_equal(reg$pval, fastH$pval, tolerance=0.01)
        fastG <- .fastG(Q, G, neig=200, ntrace=500)
        expect_equal(reg$pval, fastG$pval, tolerance=0.01)
})


test_that("fastSMMAT matches burden", {
	n <- 250

	### create a matrix of genetic variants to test.
	geno <- .testGenoMatrix(n=n, nsnp=300)
        weights <- rep(1, ncol(geno))

        ## mixed model
	nullmod <- .testNullmod(n, MM=TRUE)
        expect_message(smmat <- .testVariantSetSMMAT(nullmod, geno, weights, neig=100, verbose=TRUE), "fast_H")
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(smmat$pval_burden, burden$Score.pval)
        
        ## basic
	nullmod <- .testNullmod(n, MM=FALSE)
        expect_message(smmat <- .testVariantSetSMMAT(nullmod, geno, weights, neig=100, verbose=TRUE), "fast_H")
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(smmat$pval_burden, burden$Score.pval)
        
        ## mixed model - binary
	nullmod <- .testNullmod(n, MM=TRUE, binary=TRUE)
        expect_message(smmat <- .testVariantSetSMMAT(nullmod, geno, weights, neig=100, verbose=TRUE), "fast_H")
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(smmat$pval_burden, burden$Score.pval)
        
        ## basic - binary
	nullmod <- .testNullmod(n, MM=FALSE, binary=TRUE)
        expect_message(smmat <- .testVariantSetSMMAT(nullmod, geno, weights, neig=100, verbose=TRUE), "fast_H")
        burden <- .testVariantSetBurden(nullmod, geno, weights, burden.test="Score")
        expect_equal(smmat$pval_burden, burden$Score.pval)
})


test_that("fastSMMAT matches SMMAT", {
	n <- 250

	### create a matrix of genetic variants to test.
	geno <- .testGenoMatrix(n=n, nsnp=300)
        weights <- rep(1, ncol(geno))

        ## mixed model
	nullmod <- .testNullmod(n, MM=TRUE)
        smmat <- testVariantSet(nullmod, geno, weights, test="SMMAT")
        fastsmmat <- testVariantSet(nullmod, geno, weights, test="fastSMMAT", neig=100)
        expect_equal(smmat$pval_SMMAT, fastsmmat$pval_SMMAT, tolerance=0.01)      
})
