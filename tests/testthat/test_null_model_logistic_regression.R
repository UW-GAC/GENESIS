context("check null model logistic regression")

test_that("logistic", {
    dat <- .testNullInputs(binary=TRUE)
    
    nullmod <- .fitNullModel(dat$y, dat$X, family="binomial", verbose=FALSE)

    glm.mod <- glm(dat$y ~ -1 + dat$X, family = "binomial")

    expect_equal(nullmod$family$family, "binomial")
    expect_false(nullmod$hetResid)
    expect_equivalent(nullmod$fitted.values, fitted(glm.mod))
    expect_false(nullmod$family$mixedmodel)
    expect_equivalent(nullmod$resid.marginal, resid(glm.mod, type = "working"))
    expect_true(all(nullmod$fixef == summary(glm.mod)$coef))
    expect_true(is.null(nullmod$varComp))
    expect_null(nullmod$varCompCov)
    expect_equivalent(nullmod$betaCov, vcov(glm.mod))
    expect_equivalent(nullmod$fitted.values, fitted(glm.mod))
    expect_equal(nullmod$logLik, as.numeric(logLik(glm.mod)))
    expect_equal(nullmod$AIC, AIC(glm.mod))
    #expect_equivalent(nullmod$workingY, dat$y) this should be false
    expect_equivalent(nullmod$outcome, dat$y)
    expect_equivalent(nullmod$model.matrix, dat$X)
    expect_equivalent(diag(as.matrix(nullmod$cholSigmaInv)), sqrt(fitted(glm.mod)*(1-fitted(glm.mod))))
    expect_equal(nullmod$converged, glm.mod$converged)
    expect_null(nullmod$zeroFLAG)
    expect_equal(nullmod$RSS, 1)
    expect_true(is(nullmod, "GENESIS.nullModel"))
})
