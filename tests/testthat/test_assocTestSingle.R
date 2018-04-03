context("single variant tests")
library(SeqVarTools)

test_that("assocTestSingle", {
    svd <- .testData()
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    seqResetFilter(svd, verbose=FALSE)
    expect_equal(unique(assoc$variant.id), seqGetData(svd, "variant.id"))
    seqClose(svd)
})

test_that("assocTestSingle - sample selection", {
    svd <- .testData()
    samp <- sampleData(svd)$sample.id[sample(1:nrow(sampleData(svd)), 50)]
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), sample.id=samp, verbose=FALSE)
    expect_equal(nrow(nullmod$model.matrix), 50)
    assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    expect_equal(max(assoc$n.obs), 50)
    seqClose(svd)
})

test_that("assocTestSingle - reorder samples", {
    svd <- .testData()
    samp <- sample(sampleData(svd)$sample.id, 50)
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), sample.id=samp, verbose=FALSE)
    expect_equal(nrow(nullmod$model.matrix), 50)
    expect_equal(nullmod$sample.id, samp)
    assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    expect_equal(max(assoc$n.obs), 50)

    # check that we get same assoc results with samples in different order
    samp.sort <- sort(samp)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), sample.id=samp.sort, verbose=FALSE)
    expect_equal(nullmod$sample.id, samp.sort)
    seqResetFilter(svd, verbose=FALSE)
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    assoc2 <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    expect_equal(assoc, assoc2)
    
    seqClose(svd)
})


## test the lines of code that reorder the genotypes
test_that("reorder genotypes", {
    svd <- .testData()
    samp <- sample(sampleData(svd)$sample.id, 50)
    nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), sample.id=samp, verbose=FALSE)
    sample.index <- .setFilterNullModel(svd, nullmod, verbose=FALSE)
    geno <- expandedAltDosage(svd, use.names=TRUE, sparse=TRUE)[sample.index,,drop=FALSE]
    expect_equal(rownames(geno), samp)
    
    seqClose(svd)
})


test_that("assocTestSingle matches regression", {
    svd <- .testData()
    
    # multiallelic variants are handled differently
    snv <- isSNV(svd, biallelic=TRUE)
    seqSetFilter(svd, variant.sel=snv, verbose=FALSE)
    assoc1 <- regression(svd, outcome="outcome", covar=c("sex", "age"))
    
    nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    iterator <- SeqVarBlockIterator(svd, verbose=FALSE)
    assoc2 <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    expect_equal(nrow(assoc1), nrow(assoc2))
    expect_equal(assoc1$variant.id, assoc2$variant.id)
    expect_equal(assoc1$n, assoc2$n.obs)
    expect_equal(assoc1$freq, 1-assoc2$freq)
    ## this won't match exactly, because missing data is handled differently
    ## assocTestSingle imputes to the mean, while regression drops missing data
    #expect_equal(assoc1$Est, -assoc2$Est)
    #expect_equal(assoc1$SE, assoc2$Est.SE)
    #expect_equal(assoc1$Wald.Stat, (assoc2$Wald.Stat)^2)
    #expect_equal(assoc1$Wald.Pval, assoc2$Wald.pval, tolerance=.1)
    
    seqClose(svd)
})

test_that("assocTestSingle - GxE", {
    svd <- .testData()
    tmp <- sampleData(svd)
    tmp$env <- sample(letters[1:3], nrow(tmp), replace=TRUE)
    sampleData(svd) <- tmp
    iterator <- SeqVarBlockIterator(svd, variantBlock=2000, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age", "env"), verbose=FALSE)
    assoc <- assocTestSingle(iterator, nullmod, test="Wald", GxE="env", verbose=FALSE)
    expect_true(all(c("Est.G:envb", "SE.G:envb", "GxE.Stat") %in% names(assoc)))
    seqClose(svd)
})
