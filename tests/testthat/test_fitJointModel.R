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
  n <- 100
  nsnp <- 10
  dat <- .testNullInputs(n)
  geno <- .testGenoMatrix(n, nsnp = nsnp)
  nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose = FALSE)

  jointmod <- .fitJointModel(nullmod, geno)
  expect_is(jointmod, "list")
  expect_is(jointmod$pve, "numeric")
  expect_is(jointmod$fixef, "data.frame")
  expect_is(jointmod$covar, "matrix")
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
