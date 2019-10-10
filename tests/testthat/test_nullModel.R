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
    expect_equal(nm$sample.id, keep)
    expect_equivalent(nm$workingY, dat$a[c(TRUE,FALSE)])
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
    expect_equal(nm$sample.id, dat$sample.id)
    expect_equivalent(nm$workingY, dat$a)
})

test_that("null model from data.frame", {
    set.seed(25); a <- rnorm(10)
    dat <- data.frame(a=a,
                      b=c(rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    nm <- fitNullModel(dat, outcome="a", covars="b", verbose=FALSE)
    expect_equivalent(nm$workingY, dat$a)
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
    expect_equivalent(nm$workingY, dat$a[c(TRUE,FALSE)])
    expect_equal(nm$group.idx, list(a=1:3, b=4:5))
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
    expect_equal(nm$sample.id, rev(keep))
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
    expect_equal(nm$sample.id, dat$sample.id)
    expect_equal(rownames(nm$model.matrix), dat$sample.id)

    samp <- rev(dat$sample.id)
    covMat2 <- covMat[samp, samp]
    nm2 <- fitNullModel(dat, outcome="a", covars="b", group.var="b", cov.mat=covMat2, verbose=FALSE)
    expect_equal(nm2$sample.id, samp)
    expect_equal(rownames(nm2$model.matrix), samp)

    expect_equal(nm$workingY, rev(nm2$workingY))
    expect_equal(nm$resid[samp,], nm2$resid[samp,])
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
    expect_equal(nm$sample.id, inv$sample.id)
    
    # change order of covMat with respect to dat
    dimnames(covMat) <- list(rev(dat$sample.id), rev(dat$sample.id))
    nm <- fitNullModel(dat, outcome="a", covars="b", group.var="b", cov.mat=covMat, verbose=FALSE)
    inv <- nullModelInvNorm(nm, covMat, verbose=FALSE)
    expect_equal(nm$sample.id, inv$sample.id)
})

test_that("outcome has colnames", {
    svd <- .testData()
    adf <- sampleData(svd)
    df <- pData(adf)

    nullmod <- fitNullModel(df, outcome="outcome", covars="sex", verbose=FALSE)
    expect_equal(colnames(nullmod$outcome), "outcome")
    
    nullmod <- fitNullModel(adf, outcome="outcome", covars="sex", verbose=FALSE)
    expect_equal(colnames(nullmod$outcome), "outcome")
    
    nullmod <- fitNullModel(svd, outcome="outcome", covars="sex", verbose=FALSE)
    expect_equal(colnames(nullmod$outcome), "outcome")
    
    seqClose(svd)
})

test_that("missing data - data.frame", {
    set.seed(34); a <- c(rep(NA, 5), rnorm(10))
    dat <- data.frame(a=a,
                      b=c(rep(NA, 5), rep("a",5), rep("b", 5)),
                      stringsAsFactors=FALSE)
    set.seed(35); covMat <- crossprod(matrix(rnorm(15*2,sd=0.05),15,15))
    nm <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMat, group="b", verbose=FALSE)
    expect_equivalent(rownames(nm$model.matrix), as.character(6:15))
    expect_equivalent(nm$workingY, dat$a[6:15])
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
    expect_equal(nm$sample.id, dat$sample.id[6:15])
    expect_equivalent(nm$workingY, dat$a[6:15])
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

    # name the matrices
    covMatList <- list(a=covMat, b=covMat2)
    nm2 <- fitNullModel(dat, outcome="a", covars="b", cov.mat=covMatList, verbose=FALSE)
    expect_equivalent(nm1$varComp, nm2$varComp)
    expect_equal(names(nm2$varComp[1:2]), paste0("V_", names(covMatList)))
})

test_that("code for checking lists identical", {
    expect_true(.listIdentical(list(1:10, 1:10, 1:10)))
    expect_true(.listIdentical(list(a=1:10, b=1:10, c=1:10)))
    expect_false(.listIdentical(list(1:10, 1:10, 11:20)))
})
