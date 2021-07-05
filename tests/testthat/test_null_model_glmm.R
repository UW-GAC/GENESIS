context("check null model glmm")

test_that("glmm", {
    dat <- .testNullInputs(binary=TRUE)
    # With the test data, all variance components go to zero, so glm is used.
    # Create test data where that doesn't happen to test the glmm models.
    set.seed(1); y <- sample(c(0,1), 100, replace = T)
    set.seed(2); cov_mat <- crossprod(matrix(rnorm(100*100, sd=0.05), 100, 100))
    nullmod <- .fitNullModel(y, dat$X, covMatList=dat$cor.mat, family=binomial(), verbose=FALSE)

    expected_names <- c("model", "varComp", "varCompCov", "fixef",
                        "betaCov", "fit", "logLik", "logLikR", "niter", "AIC",
                        "model.matrix", "group.idx", "W", "cholSigmaInv", "converged",
                        "zeroFLAG", "RSS", "CX", "CXCXI", "RSS0")
    expect_true(setequal(names(nullmod), expected_names))

    # Check names of fit data frame.
    expected_names <- c("outcome", "workingY", "fitted.values", "resid.marginal",
                        "resid.conditional", "resid.PY", "resid.cholesky",
                        "linear.predictor")
    expect_true(setequal(names(nullmod$fit), expected_names))

    # Check model element
    expected_names <- c("hetResid", "family")
    expect_true(setequal(names(nullmod$model), expected_names))
    expect_false(nullmod$model$hetResid)
    expect_equal(nullmod$model$family$family, "binomial")
    expect_true(nullmod$model$family$mixedmodel)

    expect_true(is(nullmod, "GENESIS.nullMixedModel"))
    expect_true(nullmod$converged)
    expect_equivalent(nullmod$fit$outcome, y)
    expect_equivalent(nullmod$model.matrix, dat$X)
    expect_equivalent(nullmod$fit$linear.predictor, nullmod$fit$workingY - nullmod$fit$resid.conditional)
    expect_equal(nullmod$RSS, 1)
})

test_that("glmm fit as glm", {
    dat <- .testNullInputs(binary=TRUE)
    # Note that when running with test data, variance components go to zero.
    expect_message(nullmod <- .fitNullModel(dat$y, dat$X, covMatList=dat$cor.mat, family=binomial(), verbose=FALSE),
                   "using glm")

    expected_names <- c("model", "varComp", "varCompCov", "fixef",
                        "betaCov", "fit", "logLik", "AIC",
                        "model.matrix", "group.idx", "cholSigmaInv", "converged",
                        "zeroFLAG", "RSS", "CX", "CXCXI", "RSS0")
    expect_true(setequal(names(nullmod), expected_names))

    # Check names of fit data frame.
    expected_names <- c("outcome", "workingY", "fitted.values", "resid.marginal",
                        "resid.PY", "resid.cholesky")
    expect_true(setequal(names(nullmod$fit), expected_names))

    # Check model element
    expected_names <- c("hetResid", "family")
    expect_true(setequal(names(nullmod$model), expected_names))
    expect_false(nullmod$model$hetResid)
    expect_equal(nullmod$model$family$family, "binomial")
    expect_false(nullmod$model$family$mixedmodel)

    expect_true(is(nullmod, "GENESIS.nullModel"))
    expect_true(nullmod$converged)
    expect_equivalent(nullmod$fit$outcome, dat$y)
    expect_equivalent(nullmod$model.matrix, dat$X)
    expect_equal(nullmod$RSS, 1)
})
