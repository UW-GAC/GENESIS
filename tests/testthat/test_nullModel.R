context("null model")
library(Biobase)

test_that("design matrix from data.frame", {
    set.seed(20); a <- rnorm(10)
    set.seed(21); c <- sample(1:10, 10, replace=TRUE)
    dat <- data.frame(a=a,
                      b=c(rep("a",5), rep("b", 5)),
                      c=c,
                      d=rep(1, 10))
    dm <- createDesignMatrix(dat, outcome="a")
    expect_equivalent(dm$y, dat$a)
    expect_equal(ncol(dm$X), 1)
    expect_true(all(dm$X[,1] == 1))
    dm <- createDesignMatrix(dat, outcome="a", covars="b")
    expect_equal(colnames(dm$X)[-1], "bb")
    dm <- createDesignMatrix(dat, outcome="a", covars=c("b", "c", "b:c"))
    expect_equal(colnames(dm$X)[-1], c("bb", "c", "bb:c"))
    expect_message(createDesignMatrix(dat, outcome="a", covars="d"), "removed from the model")
})

test_that("design matrix with missing reference level", {
    set.seed(20); a <- rnorm(10)
    dat <- data.frame(a=a,
                      b=c("m", "n", rep("x", 4), rep("y", 4)))
    dat <- dat[3:10,]
    nm <- fitNullModel(dat, outcome="a", covars="b", verbose=FALSE)
    expect_equal(colnames(nm$model.matrix)[2], "by")
})

test_that("design matrix with collinear covariates", {
    set.seed(20); a <- rnorm(10)
    set.seed(21); c <- sample(1:10, 10, replace=TRUE)
    dat <- data.frame(a=a,
                      b=c(rep("a",5), rep("b", 5)),
                      c=c,
                      d=sample(1:10, 10, replace=TRUE),
                      stringsAsFactors = FALSE)
    dat$e <- dat$c
    dat$f <- dat$c + 2 * dat$d
    # No failure with no covariates.
     expect_silent(nm <- fitNullModel(dat, outcome="a"))
     # No failure with no colinear covariates.
     expect_silent(fitNullModel(dat, outcome="a", covars = c("b", "c", "d")))
     # Simple case where one covariate equals another.
     expect_error(fitNullModel(dat, outcome="a", covars = c("b", "c", "d", "e")),
                  "multicollinearity")
     # Slightly more complicated case involving linear combination of more than one covariate.
     expect_error(fitNullModel(dat, outcome="a", covars = c("b", "c", "d", "f")),
                  "multicollinearity")
})

test_that("null model", {
    set.seed(22); a <- rnorm(10)
    dat <- data.frame(sample.id=letters[1:10],
                      a=a,
                      b=c(rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    dat <- AnnotatedDataFrame(dat)
    keep <- dat$sample.id[c(TRUE,FALSE)]
    nm <- fitNullModel(dat, outcome="a", covars="b", sample.id=keep, verbose=FALSE)
    expect_equal(rownames(nm$model.matrix), keep)
    # Make sure that sample id was added to the fit data frame.
    expect_true("sample.id" %in% names(nm$fit))
    expect_equal(nm$fit$sample.id, keep)
    expect_equivalent(nm$fit$workingY, dat$a[c(TRUE,FALSE)])
    # Check model strings.
    expect_true("model" %in% names(nm))
    expected_names <- c("outcome", "covars", "formula", "hetResid", "family")
    expect_true(setequal(names(nm$model), expected_names))
    expect_equal(nm$model$outcome, "a")
    expect_equal(nm$model$covars, "b")
    expect_equal(nm$model$formula, "a ~ b")
    expect_false(nm$model$hetResid)

})

test_that("null model - cov.mat", {
    set.seed(23); a <- rnorm(10)
    dat <- data.frame(sample.id=letters[1:10],
                      a=a,
                      b=c(rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    dat <- AnnotatedDataFrame(dat)
    set.seed(24); covMat <- crossprod(matrix(rnorm(100,sd=0.05),10,10))
    dimnames(covMat) <- list(dat$sample.id, dat$sample.id)
    nm <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMat, verbose=FALSE)
    expect_equal(nm$fit$sample.id, dat$sample.id)
    expect_equivalent(nm$fit$workingY, dat$a)

    # Check model strings.
    expected_names <- c("outcome", "covars", "formula", "hetResid", "family")
    expect_true(setequal(names(nm$model), expected_names))
    expect_equal(nm$model$outcome, "a")
    expect_equal(nm$model$covars, "b")
    expect_equal(nm$model$formula, "a ~ b + (1|A)")
    covMatList <- list("mymat" = covMat)
    nm <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMatList, verbose=FALSE)
    expect_equal(nm$model$formula, "a ~ b + (1|mymat)")
    expect_false(nm$model$hetResid)
})

test_that("null model from data.frame", {
    set.seed(25); a <- rnorm(10)
    dat <- data.frame(a=a,
                      b=c(rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    nm <- fitNullModel(dat, outcome="a", covars="b", verbose=FALSE)
    expect_equivalent(nm$fit$workingY, dat$a)
    expect_equal(rownames(nm$model.matrix), as.character(1:nrow(dat)))
    expect_equal(rownames(nm$fit), rownames(nm$model.matrix))
    # Check model strings.
    expect_true("model" %in% names(nm))
    expected_names <- c("outcome", "covars", "formula", "hetResid", "family")
    expect_true(setequal(names(nm$model), expected_names))
    expect_equal(nm$model$outcome, "a")
    expect_equal(nm$model$covars, "b")
    expect_equal(nm$model$formula, "a ~ b")
    expect_false(nm$model$hetResid)
})

test_that("null model from data.frame with rownames", {
    set.seed(25); a <- rnorm(10)
    dat <- data.frame(a=a,
                      b=c(rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    keep <- letters[1:10]
    rownames(dat) <- keep
    nm <- fitNullModel(dat, outcome="a", covars="b", verbose=FALSE)
    expect_equivalent(nm$fit$workingY, dat$a)
    expect_equal(rownames(nm$model.matrix), keep)
    expect_equal(rownames(nm$fit), keep)
})

test_that("index list", {
    x <- rep(letters[1:3], each=3)
    expect_equal(list(a=1:3, b=4:6, c=7:9), .indexList(x))
    expect_equal(list(a=1:3), .indexList(rep("a", 3)))
})

test_that("group.var", {
    set.seed(26); a <- rnorm(10)
    dat <- data.frame(sample.id=letters[1:10],
                      a=a,
                      b=c(rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    dat <- AnnotatedDataFrame(dat)
    keep <- dat$sample.id[c(TRUE,FALSE)]
    nm <- fitNullModel(dat, outcome="a", covars="b", group.var="b", sample.id=keep, verbose=FALSE)
    expect_equal(rownames(nm$model.matrix), keep)
    expect_equivalent(nm$fit$workingY, dat$a[c(TRUE,FALSE)])
    expect_equal(nm$group.idx, list(a=1:3, b=4:5))
    # Check model strings.
    expect_true("model" %in% names(nm))
    expected_names <- c("outcome", "covars", "formula", "hetResid", "family")
    expect_true(setequal(names(nm$model), expected_names))
    expect_equal(nm$model$outcome, "a")
    expect_equal(nm$model$covars, "b")
    expect_equal(nm$model$formula, "a ~ b + var(b)")
    expect_true(nm$model$hetResid)
})

test_that("group.var is a factor", {
    set.seed(26); a <- rnorm(10)
    dat <- data.frame(sample.id=letters[1:10],
                      a=a,
                      b=factor(c(rep("a",5), rep("b", 5))),
                      stringsAsFactors=FALSE)
    dat <- AnnotatedDataFrame(dat)
    nm <- fitNullModel(dat, outcome="a", covars="b", group.var="b", verbose=FALSE)
    expect_equal(nm$group.idx, list(a=1:5, b=6:10))
    # Check model strings.
    expect_equal(nm$model$formula, "a ~ b + var(b)")

})

test_that("dimnames for cov.mat", {
    set.seed(27); a <- rnorm(10)
    dat <- data.frame(sample.id=letters[1:10],
                      a=a,
                      b=c(rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    set.seed(28); covMat <- crossprod(matrix(rnorm(100,sd=0.05),10,10))
    .checkRownames(covMat, dat)
    expect_error(.checkSampleId(covMat, dat))

    dimnames(covMat) <- list(1:10, 1:10)
    expect_error(.checkSampleId(covMat, dat))

    rownames(dat) <- dat$sample.id
    expect_error(.checkRownames(covMat, dat))

    dimnames(covMat) <- list(dat$sample.id, dat$sample.id)
    .checkSampleId(covMat, dat)
})

test_that("sample selection", {
    set.seed(28); a <- rnorm(10)
    dat <- data.frame(sample.id=letters[1:10],
                      a=a,
                      b=c(rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    dat <- AnnotatedDataFrame(dat)
    set.seed(29); covMat <- crossprod(matrix(rnorm(100,sd=0.05),10,10))
    dimnames(covMat) <- list(dat$sample.id, dat$sample.id)

    keep <- rev(dat$sample.id[c(TRUE,FALSE)])
    nm <- fitNullModel(dat, outcome="a", covars="b", group.var="b", cov.mat=covMat, sample.id=keep, verbose=FALSE)
    expect_equal(nm$fit$sample.id, rev(keep))
})

test_that("change sample order", {
    set.seed(300); a <- rnorm(10)
    dat <- data.frame(sample.id=letters[1:10],
                      a=a,
                      b=c(rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    dat <- AnnotatedDataFrame(dat)
    set.seed(301); covMat <- crossprod(matrix(rnorm(100,sd=0.05),10,10))
    dimnames(covMat) <- list(dat$sample.id, dat$sample.id)

    nm <- fitNullModel(dat, outcome="a", covars="b", group.var="b", cov.mat=covMat, verbose=FALSE)
    expect_equal(nm$fit$sample.id, dat$sample.id)
    expect_equal(rownames(nm$model.matrix), dat$sample.id)

    samp <- rev(dat$sample.id)
    covMat2 <- covMat[samp, samp]
    nm2 <- fitNullModel(dat, outcome="a", covars="b", group.var="b", cov.mat=covMat2, verbose=FALSE)
    expect_equal(nm2$fit$sample.id, samp)
    expect_equal(rownames(nm2$model.matrix), samp)

    expect_equal(nm$fit$workingY, rev(nm2$fit$workingY))
    expect_equal(nm$resid.PY[samp,], nm2$resid.PY[samp,])
    ## why are these not equal? in assocTestSingle, results are the same
    #expect_equal(nm$Ytilde[samp,], nm2$Ytilde[samp,])
    #expect_equivalent(as.matrix(nm$cholSigmaInv)[samp,samp], as.matrix(nm2$cholSigmaInv)[samp,samp])
})

test_that("inv norm", {
    set.seed(32); a <- rnorm(10)
    dat <- data.frame(sample.id=letters[1:10],
                      a=a,
                      b=c(rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    dat <- AnnotatedDataFrame(dat)
    set.seed(33); covMat <- crossprod(matrix(rnorm(100,sd=0.05),10,10, dimnames=list(1:10, 1:10)))
    dimnames(covMat) <- list(dat$sample.id, dat$sample.id)
    nm <- fitNullModel(dat, outcome="a", covars="b", group.var="b", cov.mat=covMat, verbose=FALSE)
    inv <- nullModelInvNorm(nm, covMat, verbose=FALSE)
    expect_equal(nm$fit$sample.id, inv$fit$sample.id)
    # Check model strings.
    expect_true("model" %in% names(nm))
    expected_names <- c("outcome", "covars", "formula", "hetResid", "family")
    expect_true(setequal(names(nm$model), expected_names))
    expect_equal(inv$model$outcome, "a")
    expect_equal(inv$model$covars, "b")
    expect_equal(inv$model$formula, "rankInvNorm(resid(a)) ~ b + (1|A) + var(b)")
    expect_true(inv$model$hetResid)
    expect_equal(nm$model$family, inv$model$family)

    # change order of covMat with respect to dat
    dimnames(covMat) <- list(rev(dat$sample.id), rev(dat$sample.id))
    nm <- fitNullModel(dat, outcome="a", covars="b", group.var="b", cov.mat=covMat, verbose=FALSE)
    inv <- nullModelInvNorm(nm, covMat, verbose=FALSE)
    expect_equal(nm$fit$sample.id, inv$fit$sample.id)

})

test_that("missing data - data.frame", {
    set.seed(34); a <- c(rep(NA, 5), rnorm(10))
    dat <- data.frame(a=a,
                      b=c(rep(NA, 5), rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    set.seed(35); covMat <- crossprod(matrix(rnorm(15*2,sd=0.05),15,15))
    nm <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMat, group="b", verbose=FALSE)
    expect_equivalent(rownames(nm$model.matrix), as.character(6:15))
    expect_equivalent(nm$fit$workingY, dat$a[6:15])
})

test_that("missing data - AnnotatedDataFrame", {
    set.seed(36); a <- c(rep(NA, 5), rnorm(10))
    dat <- data.frame(sample.id=letters[1:15],
                      a=a,
                      b=c(rep(NA, 5), rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    dat <- AnnotatedDataFrame(dat)
    set.seed(37); covMat <- crossprod(matrix(rnorm(15*2,sd=0.05),15,15))
    dimnames(covMat) <- list(dat$sample.id, dat$sample.id)
    nm <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMat, group="b", verbose=FALSE)
    expect_equal(nm$fit$sample.id, dat$sample.id[6:15])
    expect_equivalent(nm$fit$workingY, dat$a[6:15])
})

test_that("ScanAnnotationDataFrame", {
    set.seed(38); a <- rnorm(10)
    dat <- data.frame(sample.id=letters[1:10],
                      a=a,
                      b=c(rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    dat <- AnnotatedDataFrame(dat)
    set.seed(39); covMat <- crossprod(matrix(rnorm(100,sd=0.05),10,10))
    dimnames(covMat) <- list(dat$sample.id, dat$sample.id)
    nm1 <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMat, verbose=FALSE)

    dat <- pData(dat)
    names(dat)[1] <- "scanID"
    dat <- GWASTools::ScanAnnotationDataFrame(dat)
    nm2 <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMat, verbose=FALSE)
    expect_equal(nm1, nm2)
})

test_that("multiple matrices", {
    set.seed(40); a <- rnorm(10)
    samp <- letters[1:10]
    dat <- data.frame(sample.id=samp,
                      a=a,
                      b=c(rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    dat <- AnnotatedDataFrame(dat)
    set.seed(41); covMat <- crossprod(matrix(rnorm(15*2,sd=0.05),10,10, dimnames=list(samp,samp)))
    set.seed(42); covMat2 <- crossprod(matrix(rnorm(15*2,sd=0.05),10,10, dimnames=list(samp,samp)))
    covMatList <- list(covMat, covMat2)
    nm1 <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMatList, verbose=FALSE)
    nm <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMatList, verbose=FALSE)
    expect_equal(nm$model$formula, "a ~ b + (1|A) + (1|B)")

    # name the matrices
    covMatList <- list(m1=covMat, m2=covMat2)
    nm2 <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMatList, verbose=FALSE)
    expect_equivalent(nm1$varComp, nm2$varComp)
    expect_equal(names(nm2$varComp[1:2]), paste0("V_", names(covMatList)))
    expect_equal(nm2$model$formula, "a ~ b + (1|m1) + (1|m2)")

    # error if only one is named
    covMatList <- list(m1 = covMat, covMat2)
    expect_error(fitNullModel(dat, outcome="a", covars="b", cov.mat=covMatList, verbose=FALSE),
                 "Some names for cov.mat list are missing")

    # error if they have the same names
    covMatList <- list(m1 = covMat, m1 = covMat2)
    expect_error(fitNullModel(dat, outcome="a", covars="b", cov.mat=covMatList, verbose=FALSE),
                 "duplicated")
})

test_that("code for checking lists identical", {
    expect_true(.listIdentical(list(1:10, 1:10, 1:10)))
    expect_true(.listIdentical(list(a=1:10, b=1:10, c=1:10)))
    expect_false(.listIdentical(list(1:10, 1:10, 11:20)))
})


test_that("tibbles are supported", {
    set.seed(43); a <- rnorm(10)
    dat <- dplyr::tibble(sample.id=letters[1:10],
                         a=a,
                         b=c(rep("a",5), rep("b", 5)))
    set.seed(44); covMat <- crossprod(matrix(rnorm(100,sd=0.05),10,10))
    nm1 <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMat, verbose=FALSE)

    dat <- AnnotatedDataFrame(dat)
    dimnames(covMat) <- list(dat$sample.id, dat$sample.id)
    nm2 <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMat, verbose=FALSE)

    expect_equal(nm1$fixef, nm2$fixef)
})


test_that("data.tables are supported", {
    set.seed(43); a <- rnorm(10)
    dat <- data.table(sample.id=letters[1:10],
                      a=a,
                      b=c(rep("a",5), rep("b", 5)))
    set.seed(44); covMat <- crossprod(matrix(rnorm(100,sd=0.05),10,10))
    nm1 <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMat, verbose=FALSE)

    dat <- AnnotatedDataFrame(dat)
    dimnames(covMat) <- list(dat$sample.id, dat$sample.id)
    nm2 <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMat, verbose=FALSE)

    expect_equal(nm1$fixef, nm2$fixef)
})


test_that(".modelString", {
    outcome <- "x"
    covars <- c("a", "b")
    random <- c("r", "s")
    group_var <- "g"
    expect_equal("x ~ a + b + (1|r) + (1|s) + var(g)",
                 .modelString(outcome, covars, random, group_var))
    expect_equal("x ~ a + (1|r) + var(g)",
                 .modelString(outcome, covars[1], random[1], group_var))
    expect_equal("x ~ (1|r) + (1|s) + var(g)",
                 .modelString(outcome, NULL, random, group_var))
    expect_equal("x ~ a + b + var(g)",
                 .modelString(outcome, covars, NULL, group_var))
    expect_equal("x ~ ",
                 .modelString(outcome, NULL, NULL, NULL))

    # With inverse_normal = TRUE
    expect_equal("rankInvNorm(resid(x)) ~ a + b + (1|r) + (1|s) + var(g)",
                 .modelString(outcome, covars, random, group_var, inverse_normal = TRUE))
    expect_equal("rankInvNorm(resid(x)) ~ a + (1|r) + var(g)",
                 .modelString(outcome, covars[1], random[1], group_var, inverse_normal = TRUE))
    expect_equal("rankInvNorm(resid(x)) ~ (1|r) + (1|s) + var(g)",
                 .modelString(outcome, NULL, random, group_var, inverse_normal = TRUE))
    expect_equal("rankInvNorm(resid(x)) ~ a + b + var(g)",
                 .modelString(outcome, covars, NULL, group_var, inverse_normal = TRUE))
    expect_equal("rankInvNorm(resid(x)) ~ ",
                 .modelString(outcome, NULL, NULL, NULL, inverse_normal = TRUE))
})

test_that(".modelOutcomeString", {
  outcome <- "x"
  expect_equal(.modelOutcomeString(outcome), "x")
  expect_equal(.modelOutcomeString(outcome, inverse_normal = TRUE),  "rankInvNorm(resid(x))")
})



#if(FALSE){
test_that("update conditional model", {
  set.seed(23); a <- rnorm(100)
  set.seed(57); G <- matrix(rnorm(100, 100,1))
  dat <- data.frame(sample.id=make.unique(rep(letters, length.out = 100), sep = ""),
                    a=a,
                    b=c(rep("a",50), rep("b", 50)),
                    var_1=G[,1],
                    stringsAsFactors=FALSE)
  dat <- AnnotatedDataFrame(dat)
  set.seed(24); covMat <- crossprod(matrix(rnorm(10000,sd=0.05),100,100))
  dimnames(covMat) <- list(dat$sample.id, dat$sample.id)
  nullmod <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMat, verbose=FALSE)

  # With genotype.
  nullmod_geno <- fitNullModel(dat, outcome="a", covars=c("b", "var_1"), cov.mat=covMat, verbose=FALSE)
  # Updated null model with unnamed genotype.
  nullmod_update <- updateNullModCond(nullmod, G, covMatList=covMat, verbose=FALSE)

  expect_equivalent(nullmod_update$varComp, nullmod_geno$varComp, tolerance=1e-4)
  expect_equivalent(nullmod_update$fixef, nullmod_geno$fixef, tolerance=1e-4)
  expect_equivalent(nullmod_update$cholSigmaInv, nullmod_geno$cholSigmaInv, tolerance=1e-4)
  expect_equivalent(nullmod_update$varCompCov, nullmod_geno$varCompCov, tolerance=1e-4)
  expect_equal(nullmod_update$model$formula, nullmod_geno$model$formula)
  expect_equal(nullmod_update$model$covars, nullmod_geno$model$covars)

  # Genotype matrix with colnames.
  colnames(G) <- "variant"
  nullmod_update <- updateNullModCond(nullmod, G, covMatList=covMat, verbose=FALSE)
  expect_equal(nullmod_update$model$formula, "a ~ b + variant + (1|A)")
  expect_equal(nullmod_update$model$covars, c("b", "variant"))

})
#}

test_that(".updateNullModelFormat with current null model format", {
  set.seed(25); a <- rnorm(100)
  dat <- data.frame(sample.id=make.unique(rep(letters, length.out = 100), sep = ""),
                    a=a,
                    b=c(rep("a",50), rep("b", 50)),
                    stringsAsFactors=FALSE)
  dat <- AnnotatedDataFrame(dat)
  nullmod <- fitNullModel(dat, outcome="a", covars="b", verbose=FALSE)
  expect_silent(nm2 <- .updateNullModelFormat(nullmod))
  expect_equal(nm2, nullmod)
})

test_that(".updateNullModelFormat works with linear models", {
  nm_old <- readRDS("testdata/old_nullmod_lm.rds")
  expect_warning(nullmod <- .updateNullModelFormat(nm_old), "created with an older version")
  # Check for expected names.
  expected_names <- c("model", "varComp", "varCompCov", "fixef",
                      "betaCov", "fit", "logLik",
                      "AIC", "model.matrix",
                      "group.idx", "cholSigmaInv", "converged", "zeroFLAG",
                      "RSS", "CX", "CXCXI", "RSS0")
  expect_true(setequal(names(nullmod), expected_names))

  # Check names of fit data frame.
  expected_names <- c("outcome", "workingY", "fitted.values", "resid.marginal",
                      "resid.PY", "resid.cholesky", "sample.id")
  expect_true(setequal(names(nullmod$fit), expected_names))
  expect_equal(rownames(nullmod$fit), rownames(nullmod$model.matrix))

  # Check model element.
  expected_names <- c("hetResid", "family", "outcome")
  expect_true(setequal(names(nullmod$model), expected_names))
  expect_equal(nullmod$model$family, nm_old$family)
  expect_equal(nullmod$model$hetResid, nm_old$hetResid)
  expect_equal(nullmod$model$outcome, colnames(nm_old$outcome))
})

test_that(".updateNullModelFormat works with linear mixed models", {
  nm_old <- readRDS("testdata/old_nullmod_lmm.rds")
  expect_warning(nullmod <- .updateNullModelFormat(nm_old), "created with an older version")
  # Check for expected names.
  expected_names <- c("model", "varComp", "varCompCov", "fixef",
                      "betaCov", "fit", "logLik", "logLikR", "AIC", "model.matrix",
                      "group.idx", "W", "cholSigmaInv", "converged", "zeroFLAG",
                      "niter", "RSS", "CX", "CXCXI", "RSS0")
  expect_true(setequal(names(nullmod), expected_names))

  # Check names of fit data frame.
  expected_names <- c("outcome", "workingY", "fitted.values", "resid.marginal",
                      "resid.conditional", "resid.PY", "resid.cholesky", "sample.id",
                      "linear.predictor")
  expect_true(setequal(names(nullmod$fit), expected_names))
  expect_equal(rownames(nullmod$fit), rownames(nullmod$model.matrix))

  # Check model element.
  expected_names <- c("hetResid", "family", "outcome")
  expect_true(setequal(names(nullmod$model), expected_names))
  expect_equal(nullmod$model$family, nm_old$family)
  expect_equal(nullmod$model$hetResid, nm_old$hetResid)
  expect_equal(nullmod$model$outcome, colnames(nm_old$outcome))
})

test_that(".updateNullModelFormat works with logistic models", {
  nm_old <- readRDS("testdata/old_nullmod_glm.rds")
  expect_warning(nullmod <- .updateNullModelFormat(nm_old), "created with an older version")
  # Check for expected names.
  expected_names <- c("model", "varComp", "varCompCov", "fixef",
                      "betaCov", "fit", "logLik",
                      "AIC", "model.matrix",
                      "group.idx", "cholSigmaInv", "converged", "zeroFLAG",
                      "RSS", "CX", "CXCXI", "RSS0")
  expect_true(setequal(names(nullmod), expected_names))

  # Check names of fit data frame.
  expected_names <- c("outcome", "workingY", "fitted.values", "resid.marginal",
                      "resid.PY", "resid.cholesky", "sample.id")
  expect_true(setequal(names(nullmod$fit), expected_names))
  expect_equal(rownames(nullmod$fit), rownames(nullmod$model.matrix))

  # Check model element.
  expected_names <- c("hetResid", "family", "outcome")
  expect_true(setequal(names(nullmod$model), expected_names))
  expect_equal(nullmod$model$family, nm_old$family)
  expect_equal(nullmod$model$hetResid, nm_old$hetResid)
  expect_equal(nullmod$model$outcome, colnames(nm_old$outcome))
})

test_that(".updateNullModelFormat works with logistic mixed models", {
  nm_old <- readRDS("testdata/old_nullmod_glmm.rds")
  expect_warning(nullmod <- .updateNullModelFormat(nm_old), "created with an older version")
  # Check for expected names.
  expected_names <- c("model", "varComp", "varCompCov", "fixef",
                      "betaCov", "fit", "logLik", "logLikR", "niter", "AIC", "model.matrix",
                      "group.idx", "W", "cholSigmaInv", "converged", "zeroFLAG",
                      "RSS", "CX", "CXCXI", "RSS0")
  expect_true(setequal(names(nullmod), expected_names))

  # Check names of fit data frame.
  expected_names <- c("outcome", "workingY", "fitted.values", "resid.marginal",
                      "resid.conditional", "resid.PY", "resid.cholesky", "sample.id",
                      "linear.predictor")
  expect_true(setequal(names(nullmod$fit), expected_names))
  expect_equal(rownames(nullmod$fit), rownames(nullmod$model.matrix))

  # Check model element.
  expected_names <- c("hetResid", "family", "outcome")
  expect_true(setequal(names(nullmod$model), expected_names))
  expect_equal(nullmod$model$family, nm_old$family)
  expect_equal(nullmod$model$hetResid, nm_old$hetResid)
  expect_equal(nullmod$model$outcome, colnames(nm_old$outcome))
})

test_that(".updateNullModelFormat works with wls mixed models", {
  nm_old <- readRDS("testdata/old_nullmod_wls.rds")
  expect_warning(nullmod <- .updateNullModelFormat(nm_old), "created with an older version")
  # Check for expected names.
  expected_names <- c("model", "varComp", "varCompCov", "fixef",
                      "betaCov", "fit", "logLik", "logLikR", "niter", "AIC", "model.matrix",
                      "group.idx", "cholSigmaInv", "converged", "zeroFLAG",
                      "RSS", "CX", "CXCXI", "RSS0")
  expect_true(setequal(names(nullmod), expected_names))

  # Check names of fit data frame.
  expected_names <- c("outcome", "workingY", "fitted.values", "resid.marginal",
                      "resid.PY", "resid.cholesky", "sample.id")
  expect_true(setequal(names(nullmod$fit), expected_names))
  expect_equal(rownames(nullmod$fit), rownames(nullmod$model.matrix))

  # Check model element.
  expected_names <- c("hetResid", "family", "outcome")
  expect_true(setequal(names(nullmod$model), expected_names))
  expect_equal(nullmod$model$family, nm_old$family)
  expect_equal(nullmod$model$hetResid, nm_old$hetResid)
  expect_equal(nullmod$model$outcome, colnames(nm_old$outcome))
})

test_that(".updateNullModelFormat adds RSS0", {
  set.seed(25); a <- rnorm(100)
  dat <- data.frame(sample.id=make.unique(rep(letters, length.out = 100), sep = ""),
                    a=a,
                    b=c(rep("a",50), rep("b", 50)),
                    stringsAsFactors=FALSE)
  dat <- AnnotatedDataFrame(dat)
  nullmod <- fitNullModel(dat, outcome="a", covars="b", verbose=FALSE)
  expected_rss0 <- nullmod$RSS0
  nullmod$RSS0 <- NULL
  nm2 <- .updateNullModelFormat(nullmod)
  expect_true("RSS0" %in% names(nm2))
  expect_equal(nm2$RSS0, expected_rss0)

})

test_that("two.stage works", {
    set.seed(26); a <- rnorm(10)
    dat <- data.frame(sample.id=letters[1:10],
                      a=a,
                      b=c(rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    dat <- AnnotatedDataFrame(dat)
    nm <- fitNullModel(dat, outcome="a", covars="b", two.stage=TRUE, verbose=FALSE)
    expect_equal(nm$model$formula, "rankInvNorm(resid(a)) ~ b")
    expect_error(fitNullModel(dat, outcome="a", covars="b", two.stage=TRUE, family=binomial, verbose=FALSE), 'two stage model only applies when family is "gaussian"')
})
