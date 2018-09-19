context("check using various classes of matrices as inputs")

test_that("symmetric matrix", {
    dat <- .testNullInputs()
    cor.mat <- Matrix(dat$cor.mat, sparse=FALSE)
    
    nullmod <- .fitNullModel(dat$y, dat$X, cor.mat, verbose=FALSE)
    nullmod2 <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    expect_equal(nullmod, nullmod2)
})

test_that("packed matrix", {
    dat <- .testNullInputs()
    cor.mat <- Matrix(dat$cor.mat, sparse=FALSE)
    cor.mat <- pack(cor.mat)
    
    nullmod <- .fitNullModel(dat$y, dat$X, cor.mat, verbose=FALSE)
    nullmod2 <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    expect_equal(nullmod, nullmod2)
})

test_that("dense matrix", {
    dat <- .testNullInputs()
    cor.mat <- Matrix(dat$cor.mat, sparse=FALSE)
    cor.mat <- as(cor.mat, "dgeMatrix")
    
    nullmod <- .fitNullModel(dat$y, dat$X, cor.mat, verbose=FALSE)
    nullmod2 <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    expect_equal(nullmod, nullmod2)
})

test_that("sparse matrix", {
    dat <- .testNullInputs()
    dat$cor.mat[dat$cor.mat < 0.01] <- 0
    cor.mat <- Matrix(dat$cor.mat, sparse=TRUE)
    
    nullmod <- .fitNullModel(dat$y, dat$X, cor.mat, verbose=FALSE)
    nullmod2 <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    expect_equal(nullmod, nullmod2)
})

test_that("error on mismatched names", {
    s <- letters[1:10]
    mat1 <- matrix(1, nrow=10, ncol=10, dimnames=list(s,s))
    mat2 <- matrix(1, nrow=10, ncol=10, dimnames=list(rev(s), rev(s)))
    materr <- matrix(1, nrow=10, ncol=10, dimnames=list(s,rev(s)))
    expect_error(.covMatNames(materr))
    expect_error(.covMatNames(list(mat1, mat2)))
})

test_that("rownames only", {
    s <- letters[1:10]
    mat <- matrix(1, nrow=10, ncol=10, dimnames=list(s,NULL))
    expect_equal(s, .covMatNames(mat))
    expect_equal(s, .covMatNames(list(mat, mat)))
})

test_that("colnames only", {
    s <- letters[1:10]
    mat <- matrix(1, nrow=10, ncol=10, dimnames=list(NULL,s))
    expect_equal(s, .covMatNames(mat))
    expect_equal(s, .covMatNames(list(mat, mat)))
})
