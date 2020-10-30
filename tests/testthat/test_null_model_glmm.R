context("check null model glmm")

test_that("glmm", {
    dat <- .testNullInputs(binary=TRUE)

    nullmod <- .fitNullModel(dat$y, dat$X, covMatList=dat$cor.mat, family="binomial", verbose=FALSE)

    expected_names <- c("family", "hetResid", "varComp", "varCompCov", "fixef",
                        "betaCov", "fit", "logLik", "AIC", "model.matrix",
                        "group.idx", "cholSigmaInv", "converged", "zeroFLAG",
                        "RSS", "CX", "CXCXI", "RSS0")
    expect_true(setequal(names(nullmod), expected_names))

    # Check names of fit data frame.
    expected_names <- c("outcome", "workingY", "fitted.values", "resid.marginal",
                        "resid", "Ytilde")
    expect_true(setequal(names(nullmod$fit), expected_names))

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
    expect_equivalent(nullmod$fit$outcome, dat$y)
    expect_equivalent(nullmod$model.matrix, dat$X)
    expect_equal(nullmod$RSS, 1)
})
