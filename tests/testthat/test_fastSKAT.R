context("check fastSKAT")
library(SeqVarTools)

test_that("prep null model - no covar.matrix", {
    n <- 100
    nullmod <- .testNullmod(n, MM=FALSE)
    out <- nullModelFastSKAT(nullmod)
    expect_true(Matrix::isDiagonal(out$cholSigma))
    expect_true(is(out, "GENESIS.nullModelFastSKAT"))
    expect_true(is(out, "GENESIS.nullModel"))
    expect_equal(nullmod$model.matrix, out$model.matrix)
    expect_equal(diag(as.matrix(nullmod$cholSigmaInv)),
                 1/diag(as.matrix(out$cholSigma)))
})

test_that("prep null model - with covar.matrix", {
    n <- 100
    nullmod <- .testNullmod(n, MM=TRUE)
    out <- nullModelFastSKAT(nullmod)
    expect_true(is(out, "GENESIS.nullModelFastSKAT"))
    expect_true(is(out, "GENESIS.nullMixedModel"))
    expect_equal(nullmod$model.matrix, out$model.matrix)
    expect_equal(diag(as.matrix(nullmod$cholSigmaInv)),
                 1/diag(as.matrix(out$cholSigma)))
})

test_that("fastSKAT - rho=0", {
    n <- 100
    nullmod <- .testNullmod(n, MM=TRUE)
    out <- nullModelFastSKAT(nullmod)
    geno <- .testGenoMatrix(n)
    weights <- rep(1, ncol(geno))
    p <- .testVariantSetFastSKAT(out, geno, weights, rho=0)

    skat.p <- .testVariantSetSKAT(nullmod, geno, weights, rho=0, pval.method="davies")
    expect_equivalent(p, skat.p$pval_0, tolerance=0.01)
})

test_that("fastSKAT - rho=0.5", {
    n <- 100
    nullmod <- .testNullmod(n, MM=TRUE)
    out <- nullModelFastSKAT(nullmod)
    geno <- .testGenoMatrix(n)
    weights <- rep(1, ncol(geno))
    p <- .testVariantSetFastSKAT(out, geno, weights, rho=0.5)

    skat.p <- .testVariantSetSKAT(nullmod, geno, weights, rho=0.5, pval.method="davies")
    expect_equivalent(p, skat.p$pval_0.5, tolerance=0.01)
})

test_that("fastSKAT - rho=1", {
    n <- 100
    nullmod <- .testNullmod(n, MM=TRUE)
    out <- nullModelFastSKAT(nullmod)
    geno <- .testGenoMatrix(n)
    weights <- rep(1, ncol(geno))
    p <- .testVariantSetFastSKAT(out, geno, weights, rho=1)

    skat.p <- .testVariantSetSKAT(nullmod, geno, weights, rho=1, pval.method="davies")
    expect_equivalent(p, skat.p$pval_1, tolerance=0.01)
})

test_that("fastSKAT in assocTestAggregate", {
    svd <- .testData()
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=5e5, windowShift=2.5e5, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    
    expect_error(assocTestAggregate(iterator, nullmod, test="fastSKAT", verbose=FALSE),
                 "nullModelFastSKAT")

    nullmod <- nullModelFastSKAT(nullmod)
    assoc <- assocTestAggregate(iterator, nullmod, test="fastSKAT", verbose=FALSE)
    expect_true("pval" %in% names(assoc$results))
    
    seqClose(svd)
})
