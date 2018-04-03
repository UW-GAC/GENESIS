context("rare variant tests")
library(SeqVarTools)
library(GenomicRanges)

test_that("window", {
    svd <- .testData()
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=5e5, windowShift=2.5e5, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestAggregate(iterator, nullmod, verbose=FALSE)
    nwin <- length(variantRanges(iterator))
    expect_equal(nrow(assoc$results), nwin)
    expect_equal(length(assoc$variantInfo), nwin)
    expect_true(all(assoc$results$chr == 1))
    seqClose(svd)
})

test_that("ranges", {
    svd <- .testData()
    gr <- GRanges(seqnames=rep(1,3), ranges=IRanges(start=c(1e6, 2e6, 3e6), width=1e6))
    names(gr) <- letters[1:3]
    iterator <- SeqVarRangeIterator(svd, variantRanges=gr, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestAggregate(iterator, nullmod, verbose=FALSE)
    expect_equal(nrow(assoc$results), length(gr))
    expect_equal(length(assoc$variantInfo), length(gr))
    expect_equal(rownames(assoc$results), letters[1:3])
    expect_equal(names(assoc$variantInfo), letters[1:3])
    seqClose(svd)
})

test_that("list", {
    svd <- .testData()
    gr <- GRangesList(
        GRanges(seqnames=rep(1,2), ranges=IRanges(start=c(1e6, 3e6), width=1e6)),
        GRanges(seqnames=rep(1,2), ranges=IRanges(start=c(3e6, 34e6), width=1e6)))
    names(gr) <- letters[1:2]
    iterator <- SeqVarListIterator(svd, variantRanges=gr, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestAggregate(iterator, nullmod, verbose=FALSE)
    expect_equal(nrow(assoc$results), length(gr))
    expect_equal(length(assoc$variantInfo), length(gr))
    expect_equal(rownames(assoc$results), letters[1:2])
    expect_equal(names(assoc$variantInfo), letters[1:2])
    seqClose(svd)
})

test_that("user weights", {
    svd <- .testData()
    variant.id <- seqGetData(svd, "variant.id")
    weights <- sample(1:10, length(variant.id), replace=TRUE)
    variantData(svd) <- AnnotatedDataFrame(data.frame(variant.id, weights))
    seqSetFilterChrom(svd, include=22, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=5e5, windowShift=2.5e5, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestAggregate(iterator, nullmod, weight.user="weights", verbose=FALSE)
    tmp <- do.call(rbind, assoc$variantInfo)[,c("variant.id", "weight")]
    tmp <- tmp[!duplicated(tmp$variant.id),]
    seqSetFilterChrom(svd, include=22, verbose=FALSE)
    expect_equal(tmp$weight, variantData(svd)$weight)
    seqClose(svd)
})

test_that("exclude monomorphs", {
    svd <- .testData()
    iterator <- SeqVarBlockIterator(svd, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    freq <- seqAlleleFreq(svd)
    ind <- which(!(freq %in% c(0,1)))
    n <- length(ind)
    assoc <- assocTestAggregate(iterator, nullmod, verbose=FALSE)
    expect_equal(assoc$results$n.site, n)
    expect_equal(assoc$variantInfo[[1]]$variant.id, ind)
    seqClose(svd)
})

test_that("exclude common", {
    svd <- .testData()
    iterator <- SeqVarBlockIterator(svd, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    freq <- seqAlleleFreq(svd, ref.allele=1)
    ind <- which(!(freq %in% c(0,1)) & freq <= 0.05)
    n <- length(ind)
    assoc <- assocTestAggregate(iterator, nullmod, AF.max=0.05, verbose=FALSE)
    expect_equal(assoc$results$n.site, n)
    expect_equal(assoc$variantInfo[[1]]$variant.id, ind)
    seqClose(svd)
})

test_that("select alleles", {
    svd <- .testData()
    vi <- variantInfo(svd, expanded=TRUE)
    vi <- vi[vi$chr %in% 21:22 & vi$allele.index == 1,]
    gr <- GRanges(seqnames=vi$chr, ranges=IRanges(start=vi$pos, width=1), alt=vi$alt)
    grl <- GRangesList(gr[seqnames(gr) == 21], gr[seqnames(gr) == 22])
    iterator <- SeqVarListIterator(svd, variantRanges=grl, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestAggregate(iterator, nullmod, verbose=FALSE)
    for (i in 1:length(grl)) {
        expect_true(max(assoc$results$n.site[i]) <= length(grl[[i]]))
        expect_true(all(assoc$allele.index == 1))
    }
    seqClose(svd)
})

test_that("matchAlleles handles empty var.info", {
    svd <- .testData()
    gr <- GRanges(seqnames=1, ranges=IRanges(start=2000000, width=300000), ref="A", alt="T")
    iterator <- SeqVarRangeIterator(svd, variantRanges=gr, verbose=FALSE)
    var.info <- variantInfo(iterator)
    ind <- .matchAlleles(svd, var.info)
    expect_equal(length(ind), 0)
    seqClose(svd)
})

test_that("select alleles with empty range", {
    svd <- .testData()
    grl <- GRangesList(GRanges(seqnames=1, ranges=IRanges(start=2000000, width=1), ref="A", alt="T"))
    iterator <- SeqVarListIterator(svd, variantRanges=grl, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestAggregate(iterator, nullmod, verbose=FALSE)
    expect_equal(assoc$results$n.site, 0)
    expect_equal(length(assoc$results$variantInfo), 0)
    seqClose(svd)
})

test_that("select alleles with mismatch", {
    svd <- .testData()
    vi <- variantInfo(svd)[1:10,]
    vi$alt[2] <- vi$ref[2]
    gr <- GRanges(seqnames=vi$chr, ranges=IRanges(start=vi$pos, width=1), alt=vi$alt)
    grl <- GRangesList(gr)
    iterator <- SeqVarListIterator(svd, variantRanges=grl, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestAggregate(iterator, nullmod, verbose=FALSE)
    expect_equal(assoc$results$n.site, 9)
    expect_equal(nrow(assoc$variantInfo[[1]]), 9)
    seqClose(svd)
})
