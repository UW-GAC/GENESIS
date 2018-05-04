context("single variant tests on GenotypeIterator objects")

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


test_that("assocTestSingle - reorder samples", {
    genoData <- .testGenoData()
    samp <- as.character(sample(getScanID(genoData), 50))
    iterator <- GenotypeBlockIterator(genoData, snpBlock=1000)

    nullmod <- fitNullModel(genoData, outcome="outcome", covars="sex", sample.id=samp, verbose=FALSE)
    assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    expect_equal(nrow(nullmod$model.matrix), 50)
    expect_equal(nullmod$sample.id, samp)

    # check that we get same assoc results with samples in different order
    samp.sort <- sort(samp)
    nullmod <- fitNullModel(genoData, outcome="outcome", covars="sex", sample.id=samp.sort, verbose=FALSE)
    expect_equal(nullmod$sample.id, samp.sort)
    iterator <- GenotypeBlockIterator(genoData, snpBlock=1000)
    assoc2 <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    expect_equal(assoc, assoc2, tolerance=1e-3)
    
    close(genoData)
})
