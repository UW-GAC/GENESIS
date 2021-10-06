context("meta analysis tests")
library(SeqVarTools)
library(GenomicRanges)
library(Biobase)

BPPARAM <- BiocParallel::SerialParam()
#BPPARAM <- BiocParallel::MulticoreParam()

test_that("metaPrepScores", {
    svd <- .testData()
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=5e5, windowShift=2.5e5, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    scores <- metaPrepScores(iterator, nullmod, BPPARAM=BPPARAM, verbose=FALSE)
    nwin <- length(variantRanges(iterator))
    expect_equal(length(scores$scores.cov), nwin)
    seqClose(svd)
})
