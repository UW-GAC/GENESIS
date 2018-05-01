context("check variant set association tests")

test_that("SKAT with rho=1 matches burden", {
	n <- 100

	### create a matrix of genetic variants to test.
	geno <- .testGenoMatrix(n=n)
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


test_that("SMMAT matches burden and skat", {
	n <- 100

	### create a matrix of genetic variants to test.
	geno <- .testGenoMatrix(n=n)
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
