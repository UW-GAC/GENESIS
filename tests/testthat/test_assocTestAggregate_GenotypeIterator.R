context("aggregate variant tests on GenotypeIterator objects")
library(GWASTools)

.testSnpFilter <- function(gdsobj, breaks=100, n=10) {
    snp.index <- 1:nsnp(gdsobj)
    ind <- cut(snp.index, breaks=breaks)
    snp.list <- lapply(unique(ind), function(i) snp.index[ind == i])[1:n]
    names(snp.list) <- letters[seq_along(snp.list)]
    snp.list
}

test_that("assocTestAggregate", {
    genoData <- .testGenoData()
    covMat <- .testGenoDataGRM(genoData)
    iterator <- GenotypeIterator(genoData, snpFilter=.testSnpFilter(genoData))
    
    nullmod <- fitNullModel(genoData, outcome="outcome", covars="sex", cov.mat=covMat, verbose=FALSE)
    assoc <- assocTestAggregate(iterator, nullmod, verbose=FALSE)
    nwin <- length(snpFilter(iterator))
    expect_equal(nrow(assoc$results), nwin)
    expect_equal(length(assoc$variantInfo), nwin)
    expect_equal(rownames(assoc$results), letters[1:10])
    expect_equal(names(assoc$variantInfo), letters[1:10])
    expect_false(any(is.na(do.call(rbind, assoc$variantInfo)$variant.id)))

    close(genoData)
})

test_that("user weights", {
    genoData <- .testGenoData()
    snpID <- getSnpID(genoData)
    chromosome <- getChromosome(genoData)
    position <- getPosition(genoData)
    set.seed(8); weights <- sample(0:2, length(snpID), replace=TRUE)
    genoData@snpAnnot <- SnpAnnotationDataFrame(data.frame(snpID, chromosome, position, weights))
    iterator <- GenotypeIterator(genoData, snpFilter=.testSnpFilter(genoData))
    
    nullmod <- fitNullModel(genoData, outcome="outcome", covars="sex", verbose=FALSE)
    assoc <- assocTestAggregate(iterator, nullmod, weight.user="weights", verbose=FALSE)
    tmp <- do.call(rbind, assoc$variantInfo)
    annot <- getSnpAnnotation(genoData)
    expect_equal(tmp$weight, annot$weight[annot$snpID %in% tmp$variant.id & annot$weight > 0])

    close(genoData)
})
