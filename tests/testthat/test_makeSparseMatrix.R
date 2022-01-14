context("makeSparseMatrix tests")

test_that("no clusters", {
    x <- diag(x=1:10)
    expect_error(makeSparseMatrix(x, thresh=0.5, verbose=FALSE), "dimnames")
    colnames(x) <- 1:10
    sm <- makeSparseMatrix(x, thresh=0.5, verbose=FALSE)
    expect_true(is(sm, "sparseMatrix"))
    expect_true(isDiagonal(sm))
    expect_true(setequal(diag(x), diag(sm)))
})

test_that("two clusters", {
    x <- diag(nrow=10, ncol=10)
    rownames(x) <- 1:10
    x[2:4,2:4] <- 0.5
    x[7:9,7:9] <- 0.5
    diag(x) <- 1
    expect_message(sm <- makeSparseMatrix(x, thresh=0), "2 clusters")
    expect_true(is(sm, "sparseMatrix"))
})

test_that("sample include", {
    x <- diag(x=1:10)
    dimnames(x) <- list(1:10,1:10)
    sm <- makeSparseMatrix(x, sample.include=2:6, thresh=0, verbose=FALSE)
    expect_equal(dimnames(sm), list(as.character(2:6), as.character(2:6)))
    expect_equal(diag(sm), 2:6)
})


test_that("numeric IDs", {
    ID1 <- c(1:30,29)
    ID2 <- c(1:30,30)
    value <- c(rep(1,30),0.23)
    kin.dat <- data.frame(ID1,ID2,value)
    sm <- makeSparseMatrix(kin.dat,thresh=NULL)
    expect_true(setequal(rownames(sm), as.character(1:30)))
})
