context("single variant tests on GenotypeIterator objects")
library(GWASTools)

test_that("assocTestSingle", {
    genoData <- .testGenoData()
    covMat <- .testGenoDataGRM(genoData)
    iterator <- GenotypeBlockIterator(genoData, snpBlock=1000)
    
    nullmod <- fitNullModel(genoData, outcome="outcome", covars="sex", cov.mat=covMat, verbose=FALSE)
    assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    freq <- GWASTools::alleleFrequency(genoData)
    keep <- !is.na(freq[,"MAF"]) & freq[,"MAF"] > 0
    expect_equal(assoc$variant.id, getSnpID(genoData)[keep])

    close(genoData)
})


test_that("assocTestSingle - sample selection", {
    genoData <- .testGenoData()
    set.seed(80); samp <- getScanID(genoData)[sample(1:nscan(genoData), 50)]
    iterator <- GenotypeBlockIterator(genoData, snpBlock=1000)

    nullmod <- fitNullModel(genoData, outcome="outcome", covars="sex", sample.id=samp, verbose=FALSE)
    assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    expect_equal(max(assoc$n.obs), 50)

    close(genoData)
})


test_that("code to reorder samples works as expected", {
    genoData <- .testGenoData()
    sample.id <- getScanID(genoData)
    set.seed(81); samp.reorder <- sample(sample.id)[1:10]
    sample.index <- match(samp.reorder, sample.id)

    geno1 <- GWASTools::getGenotype(genoData, use.names=TRUE)
    geno2 <- GWASTools::getGenotypeSelection(genoData, scanID=samp.reorder, order="selection",
                                             use.names=TRUE)
    expect_equal(geno1[,sample.index], geno2)
    
    close(genoData)
})

test_that("assocTestSingle - reorder samples", {
    genoData <- .testGenoData()
    covMat <- .testGenoDataGRM(genoData)
    set.seed(82); samp <- as.character(sample(getScanID(genoData), 50))
    iterator <- GenotypeBlockIterator(genoData, snpBlock=1000)

    nullmod <- fitNullModel(genoData, outcome="outcome", cov.mat=covMat[samp,samp], verbose=FALSE)
    assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    expect_equal(nrow(nullmod$model.matrix), 50)
    expect_equal(nullmod$sample.id, samp)

    # check that we get same assoc results with samples in different order
    samp.sort <- sort(samp)
    nullmod2 <- fitNullModel(genoData, outcome="outcome", cov.mat=covMat[samp.sort,samp.sort], verbose=FALSE)
    expect_equal(nullmod2$sample.id, samp.sort)
    GWASTools::resetIterator(iterator)
    assoc2 <- assocTestSingle(iterator, nullmod2, verbose=FALSE)
    # this test may not be reliable - see test_nullModel.R
    expect_equal(assoc, assoc2)
    #expect_equal(assoc[,1:6], assoc2[,1:6])
    
    close(genoData)
})


test_that("missing sample.id in null model", {
    genoData <- .testGenoData()
    iterator <- GenotypeBlockIterator(genoData)
    samp <- pData(getScanAnnotation(genoData))[11:20,]
    nullmod <- fitNullModel(samp, outcome="outcome", covars="sex", verbose=FALSE)
    expect_false("sample.id" %in% names(nullmod))
    expect_equal(length(nullmod$outcome), 10)
    assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    expect_equal(max(assoc$n.obs), 10)
    close(genoData)
})

test_that("missing sample.id in null model - row names", {
    genoData <- .testGenoData()
    samp <- getScanAnnotation(genoData)
    sampleNames(samp) <- samp$scanID
    genoData <- GenotypeData(genoData@data, scanAnnot=samp)
    samp <- pData(getScanAnnotation(genoData))[11:20,]
    nullmod <- fitNullModel(samp, outcome="outcome", covars="sex", verbose=FALSE)
    expect_false("sample.id" %in% names(nullmod))
    expect_equal(length(nullmod$outcome), 10)
    iterator <- GenotypeBlockIterator(genoData)
    assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    expect_equal(max(assoc$n.obs), 10)
    close(genoData)
})

test_that("MatrixGenotypeReader", {
    genoData <- .testGenoData()
    geno <- GWASTools::getGenotype(genoData)
    matRdr <- MatrixGenotypeReader(geno,
                                    snpID=getSnpID(genoData),
                                    chromosome=getChromosome(genoData),
                                    position=getPosition(genoData),
                                    scanID=getScanID(genoData))
    matData <- GenotypeData(matRdr, scanAnnot=getScanAnnotation(genoData))

    iterator <- GenotypeBlockIterator(genoData, snpBlock=1000)
    nullmod <- fitNullModel(genoData, outcome="outcome", covars="sex", verbose=FALSE)
    assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)

    iterator2 <- GenotypeBlockIterator(matData, snpBlock=1000)
    nullmod <- fitNullModel(matData, outcome="outcome", covars="sex", verbose=FALSE)
    assoc2 <- assocTestSingle(iterator2, nullmod, verbose=FALSE)

    expect_equal(assoc, assoc2)

    close(genoData)
})
