context("pcrelate parallel tests")
library(SeqVarTools)

test_that("pcrelate2 - variant blocks", {
    svd <- .testData()
    mypcs <- .testPCs(svd)
    iterator <- SeqVarBlockIterator(svd, verbose=FALSE)
    myrel <- pcrelate(iterator, pcs = mypcs, num.cores=4, verbose=FALSE)
    seqResetFilter(svd, verbose=FALSE)
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    myrel2 <- pcrelate(iterator, pcs = mypcs, num.cores=4, verbose=FALSE)
    expect_equal(myrel, myrel2)
    seqResetFilter(svd, verbose=FALSE)
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    myrel3 <- pcrelate(iterator, pcs = mypcs, num.cores=1, verbose=FALSE)
    expect_equal(myrel2, myrel3)
    seqClose(svd)
})

test_that("pcrelate2 - 2 sample blocks", {
    svd <- .testData()
    mypcs <- .testPCs(svd)
    iterator <- SeqVarBlockIterator(svd, verbose=FALSE)
    myrel <- pcrelate(iterator, pcs = mypcs, small.samp.correct=FALSE, num.cores=4, verbose=FALSE)
    resetIterator(iterator, verbose=FALSE)
    myrel2 <- pcrelate(iterator, pcs = mypcs, sample.block.size=50, small.samp.correct=FALSE, num.cores=4, verbose=FALSE)
    expect_equal(myrel, myrel2)
    seqClose(svd)
})

test_that("pcrelate2 - >2 sample blocks", {
    svd <- .testData()
    mypcs <- .testPCs(svd)
    iterator <- SeqVarBlockIterator(svd, verbose=FALSE)
    myrel <- pcrelate(iterator, pcs = mypcs, small.samp.correct=FALSE, num.cores=4, verbose=FALSE)
    resetIterator(iterator, verbose=FALSE)
    myrel2 <- pcrelate(iterator, pcs = mypcs, sample.block.size=20, small.samp.correct=FALSE, num.cores=4, verbose=FALSE)
    expect_equal(myrel, myrel2)
    seqClose(svd)
})


test_that("pcrelate2 - GenotypeData - variant blocks", {
    gd <- .testGenoData()
    mypcs <- .testGenoPCs(gd)
    iterator <- GWASTools::GenotypeBlockIterator(gd)
    myrel <- pcrelate(iterator, pcs = mypcs, num.cores=4, verbose=FALSE)
    iterator <- GWASTools::GenotypeBlockIterator(gd, snpBlock=1000)
    myrel2 <- pcrelate(iterator, pcs = mypcs, num.cores=4, verbose=FALSE)
    expect_equal(myrel, myrel2)
    iterator <- GWASTools::GenotypeBlockIterator(gd, snpBlock=1000)
    myrel3 <- pcrelate(iterator, pcs = mypcs, num.cores=1, verbose=FALSE)
    expect_equal(myrel2, myrel3)
    GWASTools::close(gd)
})

test_that("pcrelate2 - GenotypeData - sample blocks", {
    gd <- .testGenoData()
    mypcs <- .testGenoPCs(gd)
    iterator <- GWASTools::GenotypeBlockIterator(gd)
    myrel <- pcrelate(iterator, pcs = mypcs, small.samp.correct=FALSE, num.cores=4, verbose=FALSE)
    iterator <- GWASTools::GenotypeBlockIterator(gd)
    myrel2 <- pcrelate(iterator, pcs = mypcs, sample.block.size=50, small.samp.correct=FALSE, num.cores=4, verbose=FALSE)
    expect_equal(myrel, myrel2)
    GWASTools::close(gd)
})

test_that("pcrelateSampBlock with 1 sample in block 1", {
    svd <- .testData()
    mypcs <- .testPCs(svd)
    iterator <- SeqVarBlockIterator(svd, verbose=FALSE)
    samp <- seqGetData(svd, "sample.id")
    beta <- calcISAFBeta(iterator, pcs=mypcs, sample.include=samp[1:10], num.cores=4, verbose=FALSE)
    resetIterator(iterator, verbose=FALSE)
    myrel <- pcrelateSampBlock(iterator, betaobj=beta, pcs=mypcs,
                               sample.include.block1=samp[1],
                               sample.include.block2=samp[2:10],
                               num.cores=4,
                               verbose=FALSE)
    expect_equal(nrow(myrel$kinBtwn), 9)
    seqClose(svd)
})
