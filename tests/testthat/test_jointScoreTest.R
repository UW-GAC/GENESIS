context("joint model")
library(Biobase)

test_that("returns expected names", {
  dat <- .testJointInputs(nsamp=100, nsnp=10)
  jointmod <- jointScoreTest(dat$nullmod, dat$geno)
  expect_is(jointmod, "list")
  expect_equal(names(jointmod), c("Joint.Stat", "Joint.pval", "Joint.PVE", "fixef", "betaCov"))
})

test_that("returns expected types", {
  dat <- .testJointInputs(nsamp=100, nsnp=10)
  jointmod <- jointScoreTest(dat$nullmod, dat$geno)
  expect_is(jointmod, "list")
  expect_is(jointmod$Joint.Stat, "numeric")
  expect_is(jointmod$Joint.pval, "numeric")
  expect_is(jointmod$Joint.PVE, "numeric")
  expect_is(jointmod$fixef, "data.frame")
  expect_is(jointmod$betaCov, "matrix")
})

test_that("works with one snp", {
  dat <- .testJointInputs(nsamp=100, nsnp=1)
  jointmod <- jointScoreTest(dat$nullmod, dat$geno)
  expect_equal(length(jointmod$Joint.Stat), 1)
  expect_equal(length(jointmod$Joint.pval), 1)
  expect_equal(length(jointmod$Joint.PVE), 1)
  expect_equal(nrow(jointmod$fixef), 1)
  expect_equal(names(jointmod$fixef), c("Est", "SE", "Stat", "pval", "PVE"))
  # Special check here - total pve should be the same as indiidual pve since its only one variant.
  expect_equal(jointmod$Joint.PVE, jointmod$fixef$PVE)
  expect_equal(nrow(jointmod$betaCov), 1)
  expect_equal(ncol(jointmod$betaCov), 1)
})

test_that("works with two snps", {
  dat <- .testJointInputs(nsamp=100, nsnp=2)
  jointmod <- jointScoreTest(dat$nullmod, dat$geno)
  expect_equal(length(jointmod$Joint.Stat), 1)
  expect_equal(length(jointmod$Joint.pval), 1)
  expect_equal(length(jointmod$Joint.PVE), 1)
  expect_equal(nrow(jointmod$fixef), 2)
  expect_equal(names(jointmod$fixef), c("Est", "SE", "Stat", "pval", "PVE"))
  expect_equal(nrow(jointmod$betaCov), 2)
  expect_equal(ncol(jointmod$betaCov), 2)
})

test_that("two snps in perfect LD", {
  dat <- .testJointInputs(nsamp=100, nsnp=2)
  dat$geno[, 2] <- dat$geno[, 1]
  expect_error(jointScoreTest(dat$nullmod, dat$geno))
})

test_that("geno matrix has no colnames", {
  dat <- .testJointInputs(nsamp=100, nsnp=2)
  jointmod <- jointScoreTest(dat$nullmod, dat$geno)
  expect_equal(rownames(jointmod$fixef), c("1", "2"))
  expect_true(is.null(rownames(jointmod$betaCov)))
  expect_true(is.null(colnames(jointmod$betaCov)))
})

test_that("uses names if geno matrix has colnames", {
  dat <- .testJointInputs(nsamp=100, nsnp=2)
  expected_names <- c("aaa", "bbb")
  colnames(dat$geno) <- expected_names
  jointmod <- jointScoreTest(dat$nullmod, dat$geno)
  expect_equal(rownames(jointmod$fixef), expected_names)
  expect_equal(rownames(jointmod$betaCov), expected_names)
  expect_equal(colnames(jointmod$betaCov), expected_names)
})

test_that("reorders samples", {
  dat <- .testJointInputs(nsamp=100, nsnp=2)
  expected_output <- jointScoreTest(dat$nullmod, dat$geno)
  idx <- sample(1:100)
  expect_equal(jointScoreTest(dat$nullmod, dat$geno[idx, ]), expected_output)
})

test_that("extra samples in geno matrix", {
  dat <- .testJointInputs(nsamp=100, nsnp=2)
  tmp <- .testGenoMatrix(n=10, nsnp=2)
  rownames(tmp) <- letters[1:10]
  extra_geno <- rbind(dat$geno, tmp)

  expected_output <- jointScoreTest(dat$nullmod, dat$geno)
  expect_equal(jointScoreTest(dat$nullmod, extra_geno), expected_output)
})

test_that("some samples missing in geno matrix", {
  dat <- .testJointInputs(nsamp=100, nsnp=2)
  tmp <- .testGenoMatrix(n=10, nsnp=2)
  idx <- sample(1:100, 10)

  expect_error(jointScoreTest(dat$nullmod, dat$geno[-idx, ]), "missing samples in genotype matrix")
})

test_that("missingness in genotype data", {
  dat <- .testJointInputs(nsamp = 100, nsnp = 2)
  dat$geno[10, 1] <- NA
  expect_error(jointScoreTest(dat$nullmod, dat$geno), "missing data")
})

test_that("works without sample.id element in null model", {
  dat <- .testJointInputs(nsamp=100, nsnp=10)
  nm <- dat$nullmod
  expected_output <- jointScoreTest(dat$nullmod, dat$geno)
  nm$fit$sample.id <- NULL
  jointmod <- jointScoreTest(nm, dat$geno)
  expect_equal(jointmod, expected_output)
})
