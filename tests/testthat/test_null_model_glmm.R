context("check null model glmm")

test_that("glmm", {
    dat <- .testNullInputs(binary=TRUE)
    
    nullmod <- .fitNullModel(dat$y, dat$X, covMatList=dat$cor.mat, family="binomial", verbose=FALSE)

    expect_equal(nullmod$family$family, "binomial")
    if (!nullmod$zeroFLAG) {
        expect_true(nullmod$family$mixedmodel)
        expect_true(is(nullmod, "GENESIS.nullMixedModel"))
    } else {
        expect_false(nullmod$family$mixedmodel)
        expect_true(is(nullmod, "GENESIS.nullModel"))
    }
    expect_false(nullmod$hetResid)
    expect_true(nullmod$converged)
    expect_equivalent(nullmod$outcome, dat$y)
    expect_equivalent(nullmod$model.matrix, dat$X)
    expect_equal(nullmod$RSS, 1)
})
