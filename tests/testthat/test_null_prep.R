context("check null model prep")

test_that("nullModelTestPrep", {
    n <- 100
    dat <- .testNullInputs(n)
    geno <- .testGenoMatrix(n)
    
    # basic
    nullmod <- .fitNullModel(dat$y, dat$X, verbose=FALSE)
    Xtilde <- calcXtilde(nullmod, geno)

    expect_equal(dim(Xtilde), c(n, ncol(geno)))
    expect_equal(dim(nullmod$Ytilde), dim(dat$y))

    # with covMatList
    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    Xtilde <- calcXtilde(nullmod, geno)

    expect_equal(dim(Xtilde), c(n, ncol(geno)))
    expect_equal(dim(nullmod$Ytilde), dim(dat$y))
})
