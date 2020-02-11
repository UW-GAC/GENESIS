context("pcair tests")

test_that("gds.class", {
    showfile.gds(closeall=TRUE)
    gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
    gds <- openfn.gds(gdsfile)
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
})


test_that("SNPGDSFileClass", {
    showfile.gds(closeall=TRUE)
    gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
    gds <- SNPRelate::snpgdsOpen(gdsfile)
    data("HapMap_ASW_MXL_KINGmat")
    mypcs <- pcair(gds, kinobj=HapMap_ASW_MXL_KINGmat, divobj=HapMap_ASW_MXL_KINGmat, verbose=FALSE)
    expect_equal(dim(mypcs$vectors), c(173,32))
    closefn.gds(gds)
})


test_that("GenotypeData", {
    showfile.gds(closeall=TRUE)
    gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
    HapMap_geno <- GWASTools::GdsGenotypeReader(gdsfile)
    data("HapMap_ASW_MXL_KINGmat")
    mypcs <- pcair(HapMap_geno, kinobj=HapMap_ASW_MXL_KINGmat, divobj=HapMap_ASW_MXL_KINGmat, verbose=FALSE)
    expect_equal(dim(mypcs$vectors), c(173,32))
    
    HapMap_genoData <- GWASTools::GenotypeData(HapMap_geno)
    mypcs <- pcair(HapMap_genoData, kinobj=HapMap_ASW_MXL_KINGmat, divobj=HapMap_ASW_MXL_KINGmat, verbose=FALSE)
    expect_equal(dim(mypcs$vectors), c(173,32))
    
    GWASTools::close(HapMap_genoData)
})


test_that("SeqVarData", {
    gds <- .testData()
    kin <- .testKing(gds)
    # weird NAs here, but we are just testing that the classes work
    mypcs <- suppressWarnings(pcair(gds, kinobj=kin, divobj=kin, verbose=FALSE))
    expect_equal(mypcs$nsamp, 90)
    seqClose(gds)
})


test_that("kinobj and divobj NULL", {
    showfile.gds(closeall=TRUE)
    gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
    gds <- openfn.gds(gdsfile)
    data("HapMap_ASW_MXL_KINGmat")

    samp <- rownames(HapMap_ASW_MXL_KINGmat)
    mypcs <- pcair(gds, kinobj = NULL, divobj = NULL, unrel.set=samp[1:100], verbose=FALSE)
    expect_equal(mypcs$rels, samp[101:173])
    expect_equal(mypcs$unrels, samp[1:100])
    
    mypcs <- pcair(gds, kinobj = NULL, divobj = NULL, unrel.set=samp[1:100], sample.include=samp[1:110], verbose=FALSE)
    expect_equal(mypcs$rels, samp[101:110])
    expect_equal(mypcs$unrels, samp[1:100])

    mypcs <- pcair(gds, kinobj = NULL, divobj = NULL, sample.include=samp[1:110], verbose=FALSE)
    expect_null(mypcs$rels)
    expect_equal(mypcs$unrels, samp[1:110])

    mypcs <- pcair(gds, verbose=FALSE)
    expect_null(mypcs$rels)
    expect_equal(mypcs$unrels, samp)
    
    closefn.gds(gds)
})


test_that("divobj NULL", {
    showfile.gds(closeall=TRUE)
    gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
    gds <- openfn.gds(gdsfile)
    data("HapMap_ASW_MXL_KINGmat")

    samp <- rownames(HapMap_ASW_MXL_KINGmat)
    mypcs <- pcair(gds, kinobj=HapMap_ASW_MXL_KINGmat, divobj = NULL, verbose=FALSE)
    expect_equal(mypcs$nsamp, 173)
    
    closefn.gds(gds)
})
