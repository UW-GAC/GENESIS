context("check null model WLS regression")

test_that("WLS", {
    dat <- .testNullInputs()

    nullmod <- .fitNullModel(dat$y, dat$X, group.idx=dat$group.idx, verbose=FALSE)

    expected_names <- c("family", "hetResid", "varComp", "varCompCov", "fixef",
                        "betaCov", "fit", "logLik", "logLikR", "AIC",
                        "model.matrix", "group.idx",
                        "cholSigmaInv", "converged", "zeroFLAG", "niter", "RSS",
                        "CX", "CXCXI", "RSS0")
    expect_true(setequal(names(nullmod), expected_names))

    # Check names of fit data frame.
    expected_names <- c("outcome", "workingY", "fitted.values", "resid.marginal",
                        "resid.conditional", "resid", "Ytilde")
    expect_true(setequal(names(nullmod$fit), expected_names))


    expect_equal(nullmod$family$family, "gaussian")
    expect_false(nullmod$family$mixedmodel)
    expect_true(nullmod$hetResid)
    expect_true(nullmod$converged)
    expect_null(nullmod$zeroFLAG)
    expect_equivalent(nullmod$fit$workingY, dat$y)
    expect_equivalent(nullmod$fit$outcome, dat$y)
    expect_equivalent(nullmod$model.matrix, dat$X)
    expect_true(is(nullmod, "GENESIS.nullModel"))
})
