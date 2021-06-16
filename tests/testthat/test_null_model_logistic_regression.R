context("check null model logistic regression")

test_that("logistic", {
    dat <- .testNullInputs(binary=TRUE)

    nullmod <- .fitNullModel(dat$y, dat$X, family=binomial(), verbose=FALSE)

    # Check for expected names.
    expected_names <- c("model", "varComp", "varCompCov", "fixef",
                        "betaCov", "fit", "logLik",
                        "AIC", "model.matrix",
                        "group.idx", "cholSigmaInv", "converged", "zeroFLAG",
                        "RSS", "CX", "CXCXI", "RSS0")
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

    glm.mod <- glm(dat$y ~ -1 + dat$X, family = binomial())

    expect_equivalent(nullmod$fit$fitted.values, fitted(glm.mod))
    expect_equivalent(nullmod$fit$resid.marginal, resid(glm.mod, type = "working"))
    expect_equivalent(as.matrix(nullmod$fixef), summary(glm.mod)$coef)
    expect_true(is.null(nullmod$varComp))
    expect_null(nullmod$varCompCov)
    expect_equivalent(nullmod$betaCov, vcov(glm.mod))
    expect_equivalent(nullmod$fit$fitted.values, fitted(glm.mod))
    expect_equal(nullmod$logLik, as.numeric(logLik(glm.mod)))
    expect_equal(nullmod$AIC, AIC(glm.mod))
    #expect_equivalent(nullmod$workingY, dat$y) this should be false
    expect_equivalent(nullmod$fit$outcome, dat$y)
    expect_equivalent(nullmod$model.matrix, dat$X)
    expect_equivalent(diag(as.matrix(nullmod$cholSigmaInv)), sqrt(fitted(glm.mod)*(1-fitted(glm.mod))))
    expect_equal(nullmod$converged, glm.mod$converged)
    expect_null(nullmod$zeroFLAG)
    expect_equal(nullmod$RSS, 1)
    expect_true(is(nullmod, "GENESIS.nullModel"))
})
