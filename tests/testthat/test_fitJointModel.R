context("joint model")
library(Biobase)

test_that("returns expected names", {
  n <- 100
  nsnp <- 10
  dat <- .testNullInputs(n)
  geno <- .testGenoMatrix(n, nsnp = nsnp)
  nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose = FALSE)

  jointmod <- .fitJointModel(nullmod, geno)
  expect_is(jointmod, "list")
  expect_equal(names(jointmod), c("pve", "fixef", "covar"))
})

test_that("returns expected types", {
})

test_that("checks input types", {
})

test_that("checks names match", {
})

test_that("works with one snp", {
})

test_that("works with two snps", {
})

test_that("two snps in perfect LD", {
})
