context("meta analysis tests")
library(SeqVarTools)
library(GenomicRanges)
library(Biobase)

BPPARAM <- BiocParallel::SerialParam()
#BPPARAM <- BiocParallel::MulticoreParam()

test_that("metaPrepScores - block", {
    svd <- .testData()
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarBlockIterator(svd, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    scores <- metaPrepScores(iterator, nullmod, BPPARAM=BPPARAM, verbose=FALSE)
    expect_true(is.data.frame(scores))
    expect_true(all(c("Score", "Score.SE") %in% names(scores)))
    seqClose(svd)
})

test_that("metaPrepScores - window", {
    svd <- .testData()
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=5e5, windowShift=2.5e5, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    scores <- metaPrepScores(iterator, nullmod, BPPARAM=BPPARAM, verbose=FALSE)
    nwin <- length(variantRanges(iterator))
    expect_equal(length(scores$scores.cov), nwin)
    seqClose(svd)
})

test_that("metaPrepScores - ranges", {
    svd <- .testData()
    gr <- GRanges(seqnames=rep(1,3), ranges=IRanges(start=c(1e6, 2e6, 3e6), width=1e6))
    names(gr) <- letters[1:3]
    iterator <- SeqVarRangeIterator(svd, variantRanges=gr, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    scores <- metaPrepScores(iterator, nullmod, BPPARAM=BPPARAM, verbose=FALSE)
    expect_equal(length(scores$scores.cov), length(gr))
    expect_equal(names(scores$scores.cov), letters[1:3])
    seqClose(svd)
})
