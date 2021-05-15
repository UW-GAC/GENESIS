context("BiocParallel tests")

test_that("assocTestSingle - SeqVarIterator", {
    svd <- .testData()
    iterator <- SeqVarTools::SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestSingle(iterator, nullmod, BPPARAM=BiocParallel::SerialParam(), verbose=FALSE)
    SeqVarTools::resetIterator(iterator, verbose=FALSE)
    assoc2 <- assocTestSingle(iterator, nullmod, BPPARAM=BiocParallel::MulticoreParam(), verbose=FALSE)
    expect_equal(assoc, assoc2)
    seqClose(svd)
})


test_that("assocTestSingle", {
    genoData <- .testGenoData()
    iterator <- GenotypeBlockIterator(genoData, snpBlock=1000)
    nullmod <- fitNullModel(genoData, outcome="outcome", covars="sex", verbose=FALSE)
    assoc <- assocTestSingle(iterator, nullmod, BPPARAM=BiocParallel::SerialParam(), verbose=FALSE)
    GWASTools::resetIterator(iterator)
    assoc2 <- assocTestSingle(iterator, nullmod, BPPARAM=BiocParallel::MulticoreParam(), verbose=FALSE)
    expect_equal(assoc, assoc2)
    close(genoData)
})
