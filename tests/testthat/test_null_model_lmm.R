context("check null model lmm")

test_that("lmm - with group", {
    dat <- .testNullInputs()
    
    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, group.idx=dat$group.idx, verbose=FALSE)

    expect_equal(nullmod$family$family, "gaussian")
    expect_true(nullmod$family$mixedmodel)
    expect_true(nullmod$hetResid)
    expect_true(nullmod$converged)
    expect_equivalent(nullmod$workingY, dat$y)
    expect_equivalent(nullmod$outcome, dat$y)
    expect_equivalent(nullmod$model.matrix, dat$X)
    expect_true(is(nullmod, "GENESIS.nullMixedModel"))
})


test_that("lmm - without group", {
    dat <- .testNullInputs()
    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)

    expect_false(nullmod$hetResid)
    expect_true(nullmod$converged)
    expect_equivalent(nullmod$workingY, dat$y)
    expect_equivalent(nullmod$outcome, dat$y)
    expect_equivalent(nullmod$model.matrix, dat$X)
    expect_true(is(nullmod, "GENESIS.nullMixedModel"))

})


if(FALSE){
test_that("update conditional model", {
    dat <- .testNullInputs()
    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, group.idx=dat$group.idx, verbose=FALSE)

    set.seed(57); G <- matrix(rnorm(100, 100,1))
    nullmod2 <- updateNullModCond(nullmod, G, covMatList=dat$cor.mat, verbose=FALSE)
    nullmod3 <- .fitNullModel(dat$y, cbind(dat$X, G), dat$cor.mat, group.idx=dat$group.idx, verbose=FALSE)

    expect_equivalent(nullmod2$varComp, nullmod3$varComp, tolerance=1e-5)
    expect_equivalent(nullmod2$fixef, nullmod3$fixef, tolerance=1e-5)
    expect_equivalent(nullmod2$cholSigmaInv, nullmod3$cholSigmaInv, tolerance=1e-5)
    expect_equivalent(nullmod2$varCompCov, nullmod3$varCompCov, tolerance=1e-5)
})
}
