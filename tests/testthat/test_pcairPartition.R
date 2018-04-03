context("pcairPartition tests")
library(GWASTools)

test_that("pcairPartition", {
    # file path to GDS file
    gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
    # read in GDS data
    HapMap_geno <- GdsGenotypeReader(filename = gdsfile)
    # create a GenotypeData class object
    HapMap_genoData <- GenotypeData(HapMap_geno)
    # load saved matrix of KING-robust estimates
    data("HapMap_ASW_MXL_KINGmat")

    mypart <- pcairPartition(kinMat = HapMap_ASW_MXL_KINGmat, divMat = HapMap_ASW_MXL_KINGmat)
    expect_equal(length(mypart$rels),76)
    expect_equal(length(mypart$unrels),97)

    newMat <- HapMap_ASW_MXL_KINGmat
    colnames(newMat) <- rownames(newMat) <- NULL
    expect_error(pcairPartition(kinMat = newMat, divMat = newMat),
                 "colnames and rownames of kinMat must be individual IDs")

    expect_error(pcairPartition(kinMat = HapMap_ASW_MXL_KINGmat, divMat = HapMap_ASW_MXL_KINGmat, unrel.set = 1:100),
                 "All of the samples in unrel.set must be in kinMat")

    close(HapMap_genoData)
})

