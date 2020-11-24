context("check null model lmm")

test_that("lmm - with group", {
    dat <- .testNullInputs()

    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, group.idx=dat$group.idx, verbose=FALSE)

    # Check for expected names.
    expected_names <- c("family", "varComp", "varCompCov", "fixef",
                        "betaCov", "fit", "logLik", "logLikR", "AIC", "model.matrix",
                        "group.idx", "cholSigmaInv", "converged", "zeroFLAG",
                        "niter", "RSS", "CX", "CXCXI", "RSS0", "model")
    expect_true(setequal(names(nullmod), expected_names))

    # Check names of fit data frame.
    expected_names <- c("outcome", "workingY", "fitted.values", "resid.marginal",
                        "resid.conditional", "resid.PY", "resid.cholesky")
    expect_true(setequal(names(nullmod$fit), expected_names))

    # Check names of model element.
    expected_names <- c("hetResid")
    expect_true(setequal(names(nullmod$model), expected_names))
    expect_true(nullmod$model$hetResid)

    expect_equal(nullmod$family$family, "gaussian")
    expect_true(nullmod$family$mixedmodel)
    expect_true(nullmod$converged)
    expect_equivalent(nullmod$fit$workingY, dat$y)
    expect_equivalent(nullmod$fit$outcome, dat$y)
    expect_equivalent(nullmod$model.matrix, dat$X)
    expect_true(is(nullmod, "GENESIS.nullMixedModel"))

})


test_that("lmm - without group", {
    dat <- .testNullInputs()
    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)

    expect_false(nullmod$model$hetResid)
    expect_true(nullmod$converged)
    expect_equivalent(nullmod$fit$workingY, dat$y)
    expect_equivalent(nullmod$fit$outcome, dat$y)
    expect_equivalent(nullmod$model.matrix, dat$X)
    expect_true(is(nullmod, "GENESIS.nullMixedModel"))

})
