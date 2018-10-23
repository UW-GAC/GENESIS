context("makeSparseMatrix tests")

test_that("no clusters", {
    x <- diag(x=1:10)
    sm <- makeSparseMatrix(x, verbose=FALSE)
    expect_true(is(sm, "sparseMatrix"))
    expect_true(isDiagonal(sm))
    expect_true(setequal(diag(x), diag(sm)))
})

test_that("two clusters", {
    x <- diag(nrow=10, ncol=10)
    x[2:4,2:4] <- 0.5
    x[7:9,7:9] <- 0.5
    diag(x) <- 1
    sm <- makeSparseMatrix(x)
    expect_true(is(sm, "sparseMatrix"))
})
