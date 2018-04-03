context("pcair tests")
library(GWASTools)

test_that("pcair", {
    # file path to GDS file
    gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
    # read in GDS data
    HapMap_geno <- GdsGenotypeReader(filename = gdsfile)
    # create a GenotypeData class object
    HapMap_genoData <- GenotypeData(HapMap_geno)
    # load saved matrix of KING-robust estimates
    data("HapMap_ASW_MXL_KINGmat")


    mypcs <- pcair(genoData = HapMap_genoData, kinMat = HapMap_ASW_MXL_KINGmat, divMat = HapMap_ASW_MXL_KINGmat)
    expect_identical(class(mypcs), "pcair")
    expect_equal(mypcs$sum.values, 100.3638, tolerance=1e-4)
    expect_equal(mypcs$kin.thresh, 2^(-11/2))
    expect_equal(mypcs$div.thresh, -2^(-11/2))
    expect_equal(mypcs$nsamp, 173)
    expect_equal(mypcs$nsnps, 19996)
    expect_equal(mypcs$MAF, 0.01)
    expect_equal(length(mypcs$rels),76)
    expect_equal(length(mypcs$unrels),97)

    expect_error(pcair(genoData = HapMap_genoData, kinMat = HapMap_ASW_MXL_KINGmat, divMat = HapMap_ASW_MXL_KINGmat, MAF=-1),
                 "MAF must be in")

    expect_error(pcair(genoData = HapMap_genoData, kinMat = HapMap_ASW_MXL_KINGmat, divMat = HapMap_ASW_MXL_KINGmat, scan.include=1:100),
                 "Not all of the scanID")

    newMat <- HapMap_ASW_MXL_KINGmat
    colnames(newMat) <- rownames(newMat) <- NULL
    expect_error(pcair(genoData = HapMap_genoData, kinMat = newMat, divMat = newMat),
                 "colnames and rownames of kinMat must be individual IDs")

    expect_error(pcair(genoData = HapMap_genoData, kinMat = HapMap_ASW_MXL_KINGmat, divMat = newMat),
                 "colnames and rownames of divMat must be individual IDs")

    expect_error(pcair(genoData = HapMap_genoData, kinMat = HapMap_ASW_MXL_KINGmat, divMat = HapMap_ASW_MXL_KINGmat, unrel.set=1:100),
                 "All of the scanIDs")

    
    close(HapMap_genoData)
})

