context("joint model")
library(Biobase)

test_that("returns expected names", {
  dat <- .testJointInputs(nsamp=100, nsnp=10)
  jointmod <- fitJointModel(dat$nullmod, dat$geno)
  expect_is(jointmod, "list")
  expect_equal(names(jointmod), c("pve", "fixef", "covar"))
})

test_that("returns expected types", {
  dat <- .testJointInputs(nsamp=100, nsnp=10)
  jointmod <- fitJointModel(dat$nullmod, dat$geno)
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
  dat <- .testJointInputs(nsamp=100, nsnp=1)
  jointmod <- fitJointModel(dat$nullmod, dat$geno)
  expect_equal(length(jointmod$pve), 1)
  expect_equal(nrow(jointmod$fixef), 1)
  expect_equal(names(jointmod$fixef), c("Est", "SE", "Stat", "pval"))
  expect_equal(nrow(jointmod$covar), 1)
  expect_equal(ncol(jointmod$covar), 1)
})

test_that("works with two snps", {
  dat <- .testJointInputs(nsamp=100, nsnp=2)
  jointmod <- fitJointModel(dat$nullmod, dat$geno)
  expect_equal(length(jointmod$pve), 1)
  expect_equal(nrow(jointmod$fixef), 2)
  expect_equal(names(jointmod$fixef), c("Est", "SE", "Stat", "pval"))
  expect_equal(nrow(jointmod$covar), 2)
  expect_equal(ncol(jointmod$covar), 2)
})

test_that("two snps in perfect LD", {
  dat <- .testJointInputs(nsamp=100, nsnp=2)
  dat$geno[, 2] <- dat$geno[, 1]
  expect_error(fitJointModel(dat$nullmod, dat$geno))
})

test_that("geno matrix has no colnames", {
  dat <- .testJointInputs(nsamp=100, nsnp=2)
  jointmod <- fitJointModel(dat$nullmod, dat$geno)
  expect_equal(rownames(jointmod$fixef), c("1", "2"))
  expect_true(is.null(rownames(jointmod$covar)))
  expect_true(is.null(colnames(jointmod$covar)))
})

test_that("uses names if geno matrix has colnames", {
  dat <- .testJointInputs(nsamp=100, nsnp=2)
  expected_names <- c("aaa", "bbb")
  colnames(dat$geno) <- expected_names
  jointmod <- fitJointModel(dat$nullmod, dat$geno)
  expect_equal(rownames(jointmod$fixef), expected_names)
  expect_equal(rownames(jointmod$covar), expected_names)
  expect_equal(colnames(jointmod$covar), expected_names)
})

test_that("reorders samples", {
  dat <- .testJointInputs(nsamp=100, nsnp=2)
  expected_output <- fitJointModel(dat$nullmod, dat$geno)
  idx <- sample(1:100)
  expect_equal(fitJointModel(dat$nullmod, dat$geno[idx, ]), expected_output)
})

test_that("extra samples in geno matrix", {
  dat <- .testJointInputs(nsamp=100, nsnp=2)
  tmp <- .testGenoMatrix(n=10, nsnp=2)
  rownames(tmp) <- letters[1:10]
  extra_geno <- rbind(dat$geno, tmp)

  expected_output <- fitJointModel(dat$nullmod, dat$geno)
  expect_equal(fitJointModel(dat$nullmod, extra_geno), expected_output)
})

test_that("some samples missing in geno matrix", {
  dat <- .testJointInputs(nsamp=100, nsnp=2)
  tmp <- .testGenoMatrix(n=10, nsnp=2)
  idx <- sample(1:100, 10)

  expect_error(fitJointModel(dat$nullmod, dat$geno[-idx, ]), "missing samples in genotype matrix")
})
