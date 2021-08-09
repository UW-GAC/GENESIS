context("BiocParallel tests")

test_that("stopOnError", {
    FUN <- function(x) {
        if (x < 5) {
            return(x)
        } else {
            stop("x too big")
        }
    }
    x <- 1:10
    i <- 1
    ITER <- function() {
        if (i <= length(x)) {
            res <- x[i]
        } else {
            res <- NULL
        }
        i <<- i + 1
        return(res)
    }
    chk <- bpiterate(ITER, FUN, BPPARAM=BiocParallel::SerialParam())
    expect_error(as.integer(unlist(chk)))
    expect_error(.stopOnError(chk), "x too big")
    
    i <- 1
    chk <- bpiterate(ITER, FUN, BPPARAM=BiocParallel::MulticoreParam())
    expect_error(.stopOnError(chk), "x too big")
    
    x <- 1:4
    i <- 1
    chk <- bpiterate(ITER, FUN, BPPARAM=BiocParallel::SerialParam())
    expect_equal(unlist(chk), 1:4)
})


test_that("assocTestSingle - SeqVarIterator", {
    svd <- .testData()
    iterator <- SeqVarTools::SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestSingle(iterator, nullmod, BPPARAM=BiocParallel::SerialParam(), verbose=FALSE)
    SeqVarTools::resetIterator(iterator, verbose=FALSE)
    assoc2 <- assocTestSingle(iterator, nullmod, BPPARAM=BiocParallel::MulticoreParam(), verbose=FALSE)
    expect_equal(assoc, assoc2)
    seqSetFilter(svd, verbose=FALSE)
    iterator <- SeqVarTools::SeqVarBlockIterator(svd, variantBlock=50, verbose=FALSE)
    assoc2 <- assocTestSingle(iterator, nullmod, BPPARAM=BiocParallel::MulticoreParam(), verbose=FALSE)
    expect_equal(assoc, assoc2)
    seqClose(svd)
})


test_that("assocTestSingle - GenotypeIterator", {
    genoData <- .testGenoData()
    iterator <- GWASTools::GenotypeBlockIterator(genoData, snpBlock=1000)
    nullmod <- fitNullModel(genoData, outcome="outcome", covars="sex", verbose=FALSE)
    assoc <- assocTestSingle(iterator, nullmod, BPPARAM=BiocParallel::SerialParam(), verbose=FALSE)
    GWASTools::resetIterator(iterator)
    assoc2 <- assocTestSingle(iterator, nullmod, BPPARAM=BiocParallel::MulticoreParam(), verbose=FALSE)
    expect_equal(assoc, assoc2)
    close(genoData)
})

test_that("assocTestAggreagate - SeqVarIterator", {
    svd <- .testData()
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=5e5, windowShift=2.5e5, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestAggregate(iterator, nullmod, BPPARAM=BiocParallel::SerialParam(), verbose=FALSE)
    SeqVarTools::resetIterator(iterator, verbose=FALSE)
    assoc2 <- assocTestAggregate(iterator, nullmod, BPPARAM=BiocParallel::MulticoreParam(), verbose=FALSE)
    expect_equal(assoc, assoc2)
    seqClose(svd)
})

.testSnpFilter <- function(gdsobj, breaks=100, n=10) {
    snp.index <- 1:GWASTools::nsnp(gdsobj)
    ind <- cut(snp.index, breaks=breaks)
    snp.list <- lapply(unique(ind), function(i) snp.index[ind == i])[1:n]
    names(snp.list) <- letters[seq_along(snp.list)]
    snp.list
}

test_that("assocTestAggregate - GenotypeIterator", {
    genoData <- .testGenoData()
    iterator <- GWASTools::GenotypeIterator(genoData, snpFilter=.testSnpFilter(genoData))
    nullmod <- fitNullModel(genoData, outcome="outcome", covars="sex", verbose=FALSE)
    assoc <- assocTestAggregate(iterator, nullmod, BPPARAM=BiocParallel::SerialParam(), verbose=FALSE)
    GWASTools::resetIterator(iterator)
    assoc2 <- assocTestAggregate(iterator, nullmod, BPPARAM=BiocParallel::MulticoreParam(), verbose=FALSE)
    expect_equal(assoc, assoc2)
    close(genoData)
})
