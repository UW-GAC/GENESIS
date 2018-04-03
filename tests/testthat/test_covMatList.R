context("check using various classes of matrices as inputs")

test_that("symmetric matrix", {
    dat <- .testNullInputs()
    cor.mat <- Matrix(dat$cor.mat, sparse=FALSE)
    
    nullmod <- fitNullMod(dat$y, dat$X, cor.mat, verbose=FALSE)
    nullmod2 <- fitNullMod(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    expect_equal(nullmod, nullmod2)
})

test_that("packed matrix", {
    dat <- .testNullInputs()
    cor.mat <- Matrix(dat$cor.mat, sparse=FALSE)
    cor.mat <- pack(cor.mat)
    
    nullmod <- fitNullMod(dat$y, dat$X, cor.mat, verbose=FALSE)
    nullmod2 <- fitNullMod(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    expect_equal(nullmod, nullmod2)
})

test_that("dense matrix", {
    dat <- .testNullInputs()
    cor.mat <- Matrix(dat$cor.mat, sparse=FALSE)
    cor.mat <- as(cor.mat, "dgeMatrix")
    
    nullmod <- fitNullMod(dat$y, dat$X, cor.mat, verbose=FALSE)
    nullmod2 <- fitNullMod(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    expect_equal(nullmod, nullmod2)
})

test_that("sparse matrix", {
    dat <- .testNullInputs()
    dat$cor.mat[dat$cor.mat < 0.01] <- 0
    cor.mat <- Matrix(dat$cor.mat, sparse=TRUE)
    
    nullmod <- fitNullMod(dat$y, dat$X, cor.mat, verbose=FALSE)
    nullmod2 <- fitNullMod(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    expect_equal(nullmod, nullmod2)
})
