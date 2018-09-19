context("check null model WLS regression")

test_that("WLS", {
    dat <- .testNullInputs()

    nullmod <- .fitNullModel(dat$y, dat$X, group.idx=dat$group.idx, verbose=FALSE)

    expect_equal(nullmod$family$family, "gaussian")
    expect_false(nullmod$family$mixedmodel)
    expect_true(nullmod$hetResid)
    expect_true(nullmod$converged)
    expect_null(nullmod$zeroFLAG) 
    expect_equivalent(nullmod$workingY, dat$y)
    expect_equivalent(nullmod$outcome, dat$y)
    expect_equivalent(nullmod$model.matrix, dat$X)
    expect_true(is(nullmod, "GENESIS.nullModel"))
})
