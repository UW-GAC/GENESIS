context("check that null model has expected matrix or Matrix objects")

test_that("gaussian, no covMat, no group", {
    dat <- .testNullInputs()

    nullmod <- .fitNullModel(dat$y, dat$X, verbose=FALSE)
    expect_true(is.numeric(nullmod$cholSigmaInv))
    expect_equal(length(nullmod$cholSigmaInv), 1)
    #expect_true(is.matrix(nullmod$Ytilde))
    expect_true(is.matrix(nullmod$CX))
    expect_true(is.matrix(nullmod$CXCXI))
    #expect_true(is.numeric(nullmod$resid))
    expect_true(is.matrix(nullmod$model.matrix))
})

test_that("gaussian, no covMat, group", {
    dat <- .testNullInputs()

    nullmod <- .fitNullModel(dat$y, dat$X, group.idx=dat$group.idx, verbose=FALSE)
    expect_true(is(nullmod$cholSigmaInv, "Matrix"))
    expect_true(isDiagonal(nullmod$cholSigmaInv))
    #expect_true(is(nullmod$Ytilde, "Matrix"))
    expect_true(is(nullmod$CX, "Matrix"))
    expect_true(is(nullmod$CXCXI, "Matrix"))
    #expect_true(is(nullmod$resid, "Matrix"))
    expect_true(is.matrix(nullmod$model.matrix))
})

test_that("gaussian, covMat - matrix, no group", {
    dat <- .testNullInputs()

    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    expect_true(is.matrix(nullmod$cholSigmaInv))
    #expect_true(is.matrix(nullmod$Ytilde))
    expect_true(is.matrix(nullmod$CX))
    expect_true(is.matrix(nullmod$CXCXI))
    #expect_true(is.numeric(nullmod$resid))
    expect_true(is.matrix(nullmod$model.matrix))
})

test_that("gaussian, covMat - Matrix, no group", {
    dat <- .testNullInputs()

    nullmod <- .fitNullModel(dat$y, dat$X, Matrix(dat$cor.mat), verbose=FALSE)
    expect_true(is(nullmod$cholSigmaInv, "Matrix"))
    #expect_true(is(nullmod$Ytilde, "Matrix"))
    expect_true(is(nullmod$CX, "Matrix"))
    expect_true(is(nullmod$CXCXI, "Matrix"))
    #expect_true(is(nullmod$resid, "Matrix"))
    expect_true(is.matrix(nullmod$model.matrix))
})

test_that("gaussian, covMat - matrix, group", {
    dat <- .testNullInputs()

    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, group.idx=dat$group.idx, verbose=FALSE)
    expect_true(is.matrix(nullmod$cholSigmaInv))
    #expect_true(is.matrix(nullmod$Ytilde))
    expect_true(is.matrix(nullmod$CX))
    expect_true(is.matrix(nullmod$CXCXI))
    #expect_true(is.numeric(nullmod$resid))
    expect_true(is.matrix(nullmod$model.matrix))
})

test_that("gaussian, covMat - Matrix, group", {
    dat <- .testNullInputs()

    nullmod <- .fitNullModel(dat$y, dat$X, Matrix(dat$cor.mat), group.idx=dat$group.idx, verbose=FALSE)
    expect_true(is(nullmod$cholSigmaInv, "Matrix"))
    #expect_true(is(nullmod$Ytilde, "Matrix"))
    expect_true(is(nullmod$CX, "Matrix"))
    expect_true(is(nullmod$CXCXI, "Matrix"))
    #expect_true(is(nullmod$resid, "Matrix"))
    expect_true(is.matrix(nullmod$model.matrix))
})

test_that("binary, no covMat", {
    dat <- .testNullInputs(binary=TRUE)

    nullmod <- .fitNullModel(dat$y, dat$X, family=binomial(), verbose=FALSE)
    expect_true(isDiagonal(nullmod$cholSigmaInv))
    #expect_true(is(nullmod$Ytilde, "Matrix"))
    expect_true(is(nullmod$CX, "Matrix"))
    expect_true(is(nullmod$CXCXI, "Matrix"))
    #expect_true(is(nullmod$resid, "Matrix"))
    expect_true(is.matrix(nullmod$model.matrix))
})

test_that("binary, covMat - matrix", {
    dat <- .testNullInputs(binary=TRUE)

    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, family=binomial(), verbose=FALSE)
    if (nullmod$zeroFLAG) {
        expect_true(isDiagonal(nullmod$cholSigmaInv))
        #expect_true(is(nullmod$Ytilde, "Matrix"))
        expect_true(is(nullmod$CX, "Matrix"))
        expect_true(is(nullmod$CXCXI, "Matrix"))
        #expect_true(is(nullmod$resid, "Matrix"))
    } else {
        expect_true(is.matrix(nullmod$cholSigmaInv))
        #expect_true(is.matrix(nullmod$Ytilde))
        expect_true(is.matrix(nullmod$CX))
        expect_true(is.matrix(nullmod$CXCXI))
        #expect_true(is.numeric(nullmod$resid))
    }
    expect_true(is.matrix(nullmod$model.matrix))
})

test_that("binary, covMat - Matrix", {
    dat <- .testNullInputs(binary=TRUE)

    nullmod <- .fitNullModel(dat$y, dat$X, Matrix(dat$cor.mat), family=binomial(), verbose=FALSE)
    expect_true(is(nullmod$cholSigmaInv, "Matrix"))
    #expect_true(is(nullmod$Ytilde, "Matrix"))
    expect_true(is(nullmod$CX, "Matrix"))
    expect_true(is(nullmod$CXCXI, "Matrix"))
    #expect_true(is(nullmod$resid, "Matrix"))
    if (nullmod$zeroFLAG) {
        expect_true(isDiagonal(nullmod$cholSigmaInv))
    }
    expect_true(is.matrix(nullmod$model.matrix))
})

test_that("convert to Matrix", {
    dat <- .testNullInputs()

    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    expect_true(is.matrix(nullmod$cholSigmaInv))

    nullmod <- .nullModelAsMatrix(nullmod)
    expect_true(is(nullmod$cholSigmaInv, "Matrix"))
    #expect_true(is(nullmod$Ytilde, "Matrix"))
    expect_true(is(nullmod$CX, "Matrix"))
    expect_true(is(nullmod$CXCXI, "Matrix"))
    #expect_true(is(nullmod$resid, "Matrix"))
})
