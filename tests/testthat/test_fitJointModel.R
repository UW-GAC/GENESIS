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
  # test in wrapper function
})

test_that("checks names match", {
  # test in wrapper function
})

test_that("works with one snp", {
  n <- 100
  nsnp <- 1
  dat <- .testNullInputs(n)
  geno <- .testGenoMatrix(n, nsnp = nsnp)
  nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose = FALSE)

  jointmod <- .fitJointModel(nullmod, geno)
  expect_equal(length(jointmod$pve), 1)
  expect_equal(nrow(jointmod$fixef), nsnp)
  expect_equal(names(jointmod$fixef), c("Est", "SE", "Stat", "pval"))
  expect_equal(nrow(jointmod$covar), nsnp)
  expect_equal(ncol(jointmod$covar), nsnp)
})

test_that("works with two snps", {
  n <- 100
  nsnp <- 2
  dat <- .testNullInputs(n)
  geno <- .testGenoMatrix(n, nsnp = nsnp)
  nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose = FALSE)

  jointmod <- .fitJointModel(nullmod, geno)
  expect_equal(length(jointmod$pve), 1)
  expect_equal(nrow(jointmod$fixef), nsnp)
  expect_equal(names(jointmod$fixef), c("Est", "SE", "Stat", "pval"))
  expect_equal(nrow(jointmod$covar), nsnp)
  expect_equal(ncol(jointmod$covar), nsnp)
  fail('add checks on names')
})

test_that("two snps in perfect LD", {
  n <- 100
  nsnp <- 2
  dat <- .testNullInputs(n)
  geno <- .testGenoMatrix(n, nsnp = nsnp)
  geno[, 2] <- geno[, 1]
  nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose = FALSE)

  expect_error(.fitJointModel(nullmod, geno))
})
