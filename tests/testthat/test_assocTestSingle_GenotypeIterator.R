context("single variant tests on GenotypeIterator objects")
library(GWASTools)

test_that("assocTestSingle", {
    genoData <- .testGenoData()
    covMat <- .testGenoDataGRM(genoData)
    iterator <- GenotypeBlockIterator(genoData, snpBlock=1000)
    
    #nullmod <- fitNullModel(scanAnnot, outcome="status", covars="sex", cov.mat=covMat, family="binomial", verbose=FALSE)
    nullmod <- fitNullModel(genoData, outcome="outcome", covars="sex", cov.mat=covMat, verbose=FALSE)
    assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    expect_equal(assoc$variant.id, getSnpID(genoData))

    close(genoData)
})


test_that("assocTestSingle - sample selection", {
    genoData <- .testGenoData()
    samp <- getScanID(genoData)[sample(1:nscan(genoData), 50)]
    iterator <- GenotypeBlockIterator(genoData, snpBlock=1000)

    nullmod <- fitNullModel(genoData, outcome="outcome", covars="sex", sample.id=samp, verbose=FALSE)
    assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    expect_equal(max(assoc$n.obs), 50)

    close(genoData)
})


test_that("code to reorder samples works as expected", {
    genoData <- .testGenoData()
    sample.id <- getScanID(genoData)
    samp.reorder <- sample(sample.id)[1:10]
    sample.index <- match(samp.reorder, sample.id)

    geno1 <- GWASTools::getGenotype(genoData, use.names=TRUE)
    geno2 <- GWASTools::getGenotypeSelection(genoData, scanID=samp.reorder, order="selection",
                                             use.names=TRUE)
    expect_equal(geno1[,sample.index], geno2)
    
    close(genoData)
})

test_that("assocTestSingle - reorder samples", {
    genoData <- .testGenoData()
    samp <- as.character(sample(getScanID(genoData), 50))
    iterator <- GenotypeBlockIterator(genoData, snpBlock=1000)

    nullmod <- fitNullModel(genoData, outcome="outcome", sample.id=samp, verbose=FALSE)
    assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    expect_equal(nrow(nullmod$model.matrix), 50)
    expect_equal(nullmod$sample.id, samp)

    # check that we get same assoc results with samples in different order
    samp.sort <- sort(samp)
    nullmod <- fitNullModel(genoData, outcome="outcome", sample.id=samp.sort, verbose=FALSE)
    expect_equal(nullmod$sample.id, samp.sort)
    iterator <- GenotypeBlockIterator(genoData, snpBlock=1000)
    assoc2 <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    expect_equal(assoc, assoc2, tolerance=1e-3)
    
    close(genoData)
})
