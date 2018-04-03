context("check null model linear regression")

test_that("linear regression", {
    n <- 100
    dat <- .testNullInputs(n)

    nullmod <- fitNullMod(dat$y, dat$X, verbose=FALSE)

    lm.mod <- lm(dat$y ~ -1 + dat$X)

    expect_equal(nullmod$fitted.values, fitted(lm.mod))
    expect_equal(nullmod$family$family, "gaussian")
    expect_false(nullmod$family$mixedmodel)
    expect_false(nullmod$hetResid)
    expect_equivalent(nullmod$resid.marginal, lm.mod$resid)
    expect_true(all(nullmod$fixef == summary(lm.mod)$coef))
    expect_equal(nullmod$varComp, summary(lm.mod)$sigma^2)
    expect_null(nullmod$varCompCov)
    expect_equivalent(nullmod$betaCov, vcov(lm.mod))
    expect_equivalent(nullmod$fitted.values, fitted(lm.mod))
    expect_equal(nullmod$logLik, as.numeric(logLik(lm.mod)))
    expect_equal(nullmod$AIC, AIC(lm.mod))
    expect_equivalent(nullmod$workingY, dat$y)
    expect_equivalent(nullmod$outcome, dat$y)
    expect_equivalent(nullmod$model.matrix, dat$X)
    expect_equal(nullmod$cholSigmaInv, 1/summary(lm.mod)$sigma)
    expect_true(nullmod$converged)
    expect_null(nullmod$zeroFLAG)
    expect_equal(nullmod$RSS, sum(lm.mod$resid^2)/(summary(lm.mod)$sigma^2*(n - ncol(dat$X))))
})

