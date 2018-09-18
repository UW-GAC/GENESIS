context("pcair tests")
#library(GWASTools)

test_that("pcair", {
    # file path to GDS file
    gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
    gds <- SNPRelate::snpgdsOpen(gdsfile)
    # read in GDS data
    #HapMap_geno <- GdsGenotypeReader(filename = gdsfile)
    # create a GenotypeData class object
    #HapMap_genoData <- GenotypeData(HapMap_geno)
    # load saved matrix of KING-robust estimates
    data("HapMap_ASW_MXL_KINGmat")

    mypcs <- pcair(gds, kinobj=HapMap_ASW_MXL_KINGmat, divobj=HapMap_ASW_MXL_KINGmat, eigen.cnt=20, verbose=FALSE)
    expect_identical(class(mypcs), "pcair")
    expect_equal(dim(mypcs$vectors), c(173,20))
    expect_equal(length(mypcs$values), 20)
    expect_equal(length(mypcs$rels),76)
    expect_equal(length(mypcs$unrels),97)
    expect_equal(mypcs$kin.thresh, 2^(-11/2))
    expect_equal(mypcs$div.thresh, -2^(-11/2))
    expect_equal(mypcs$nsamp, 173)

    closefn.gds(gds)
    #close(HapMap_genoData)
})

