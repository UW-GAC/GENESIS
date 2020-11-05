context("check null model prep")

test_that("nullModelTestPrep", {
    n <- 100
    dat <- .testNullInputs(n)
    geno <- .testGenoMatrix(n)



    # basic
    nullmod <- .fitNullModel(dat$y, dat$X, verbose=FALSE)
    Xtilde <- calcGtilde(nullmod, geno)

    expect_equal(dim(Xtilde), c(n, ncol(geno)))
    expect_equal(length(nullmod$fit$resid.cholesky), length(dat$y))

    # Adds expected columns to fit data frame.
    expect_true("resid" %in% names(nullmod$fit))
    expect_true("resid.cholesky" %in% names(nullmod$fit))

    # with covMatList
    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    Xtilde <- calcGtilde(nullmod, geno)

    expect_equal(dim(Xtilde), c(n, ncol(geno)))
    expect_equal(length(nullmod$fit$resid.cholesky), length(dat$y))
})
