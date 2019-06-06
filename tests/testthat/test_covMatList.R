context("check using various classes of matrices as inputs")

.compareMats <- function(nm1, nm2) {
    for (m in c("model.matrix", "cholSigmaInv", "Ytilde", "resid", "CX", "CXCXI")) {
        expect_equivalent(as.matrix(nm1[[m]]), as.matrix(nm2[[m]]))
    }
}

test_that("symmetric matrix", {
    dat <- .testNullInputs()
    cor.mat <- Matrix(dat$cor.mat, sparse=FALSE)
    
    nullmod <- .fitNullModel(dat$y, dat$X, cor.mat, verbose=FALSE)
    nullmod2 <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    .compareMats(nullmod, nullmod2)
})

test_that("packed matrix", {
    dat <- .testNullInputs()
    cor.mat <- Matrix(dat$cor.mat, sparse=FALSE)
    cor.mat <- pack(cor.mat)
    
    nullmod <- .fitNullModel(dat$y, dat$X, cor.mat, verbose=FALSE)
    nullmod2 <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    .compareMats(nullmod, nullmod2)
})

test_that("dense matrix", {
    dat <- .testNullInputs()
    cor.mat <- Matrix(dat$cor.mat, sparse=FALSE)
    cor.mat <- as(cor.mat, "dgeMatrix")
    
    nullmod <- .fitNullModel(dat$y, dat$X, cor.mat, verbose=FALSE)
    nullmod2 <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    .compareMats(nullmod, nullmod2)
})

test_that("sparse matrix", {
    dat <- .testNullInputs()
    dat$cor.mat[dat$cor.mat < 0.01] <- 0
    cor.mat <- Matrix(dat$cor.mat, sparse=TRUE)
    
    nullmod <- .fitNullModel(dat$y, dat$X, cor.mat, verbose=FALSE)
    nullmod2 <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    .compareMats(nullmod, nullmod2)
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

test_that("multiple matrices", {
    dat <- .testNullInputs()
    mat1 <- dat$cor.mat
    n <- length(dat$y)
    set.seed(1); sqrt.cor.mat <- matrix(rnorm(n*n, sd = 0.05),n,n, dimnames=list(1:n, 1:n))
    mat2 <- crossprod(sqrt.cor.mat)
    covMatList <- list(mat1, mat2)
    expect_true(all(sapply(.checkMatrixType(covMatList), is.matrix)))
    nullmod <- .fitNullModel(dat$y, dat$X, covMatList, verbose=FALSE)
    expect_true(is.matrix(nullmod$cholSigmaInv))
    
    covMatList <- list(mat1, Matrix(mat2))
    expect_true(all(sapply(.checkMatrixType(covMatList), is, "Matrix")))
    nullmod <- .fitNullModel(dat$y, dat$X, covMatList, verbose=FALSE)
    expect_true(is(nullmod$cholSigmaInv, "Matrix"))
    
    covMatList <- list(Matrix(mat1), Matrix(mat2))
    expect_true(all(sapply(.checkMatrixType(covMatList), is, "Matrix")))
    nullmod <- .fitNullModel(dat$y, dat$X, covMatList, verbose=FALSE)
    expect_true(is(nullmod$cholSigmaInv, "Matrix"))
})
