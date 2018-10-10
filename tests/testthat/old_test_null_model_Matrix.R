context("check that null model has Matrix objects")

test_that("gaussian, no covMat, no group", {
    dat <- .testNullInputs()

    nullmod <- .fitNullModel(dat$y, dat$X, verbose=FALSE)
    expect_true(is.numeric(nullmod$cholSigmaInv))
    expect_equal(length(nullmod$cholSigmaInv), 1)
    expect_equal(attr(class(nullmod$Ytilde), "package"), "Matrix")
    expect_equal(attr(class(nullmod$CX), "package"), "Matrix")
    expect_equal(attr(class(nullmod$CXCXI), "package"), "Matrix")
    expect_true(is.numeric(nullmod$resid))
    expect_true(is.matrix(nullmod$model.matrix))
})

    
test_that("gaussian, no covMat, group", {
    dat <- .testNullInputs()

    nullmod <- .fitNullModel(dat$y, dat$X, group.idx=dat$group.idx, verbose=FALSE)
    expect_equal(attr(class(nullmod$cholSigmaInv), "package"), "Matrix")
    expect_true(isDiagonal(nullmod$cholSigmaInv))
    expect_equal(attr(class(nullmod$Ytilde), "package"), "Matrix")
    expect_equal(attr(class(nullmod$CX), "package"), "Matrix")
    expect_equal(attr(class(nullmod$CXCXI), "package"), "Matrix")
    expect_equal(attr(class(nullmod$resid), "package"), "Matrix")
    expect_true(is.matrix(nullmod$model.matrix))
})

test_that("gaussian, covMat, no group", {
    dat <- .testNullInputs()

    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    expect_equal(attr(class(nullmod$cholSigmaInv), "package"), "Matrix")
    expect_equal(attr(class(nullmod$Ytilde), "package"), "Matrix")
    expect_equal(attr(class(nullmod$CX), "package"), "Matrix")
    expect_equal(attr(class(nullmod$CXCXI), "package"), "Matrix")
    expect_equal(attr(class(nullmod$resid), "package"), "Matrix")
    expect_true(is.matrix(nullmod$model.matrix))
})
    
test_that("gaussian, covMat, group", {
    dat <- .testNullInputs()

    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, group.idx=dat$group.idx, verbose=FALSE)
    expect_equal(attr(class(nullmod$cholSigmaInv), "package"), "Matrix")
    expect_equal(attr(class(nullmod$Ytilde), "package"), "Matrix")
    expect_equal(attr(class(nullmod$CX), "package"), "Matrix")
    expect_equal(attr(class(nullmod$CXCXI), "package"), "Matrix")
    expect_equal(attr(class(nullmod$resid), "package"), "Matrix")
    expect_true(is.matrix(nullmod$model.matrix))
})

test_that("binary, no covMat", {
    dat <- .testNullInputs(binary=TRUE)

    nullmod <- .fitNullModel(dat$y, dat$X, family="binomial", verbose=FALSE)
    expect_equal(attr(class(nullmod$cholSigmaInv), "package"), "Matrix")
    expect_true(isDiagonal(nullmod$cholSigmaInv))
    expect_equal(attr(class(nullmod$Ytilde), "package"), "Matrix")
    expect_equal(attr(class(nullmod$CX), "package"), "Matrix")
    expect_equal(attr(class(nullmod$CXCXI), "package"), "Matrix")
    expect_true(is.numeric(nullmod$resid))
    expect_true(is.matrix(nullmod$model.matrix))
})

test_that("binary, covMat", {
    dat <- .testNullInputs(binary=TRUE)

    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, family="binomial", verbose=FALSE)
    expect_equal(attr(class(nullmod$cholSigmaInv), "package"), "Matrix")
    expect_equal(attr(class(nullmod$Ytilde), "package"), "Matrix")
    expect_equal(attr(class(nullmod$CX), "package"), "Matrix")
    expect_equal(attr(class(nullmod$CXCXI), "package"), "Matrix")
    if (!nullmod$zeroFLAG) expect_equal(attr(class(nullmod$resid), "package"), "Matrix")
    if (nullmod$zeroFLAG) {
        expect_true(isDiagonal(nullmod$cholSigmaInv))
        expect_true(is.numeric(nullmod$resid))
    }   
    expect_true(is.matrix(nullmod$model.matrix))
})
