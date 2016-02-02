test_pcair <- function(){
    # file path to GDS file
    gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
    # read in GDS data
    HapMap_geno <- GdsGenotypeReader(filename = gdsfile)
    # create a GenotypeData class object
    HapMap_genoData <- GenotypeData(HapMap_geno)
    # load saved matrix of KING-robust estimates
    data("HapMap_ASW_MXL_KINGmat")


    mypcs <- pcair(genoData = HapMap_genoData, kinMat = HapMap_ASW_MXL_KINGmat, divMat = HapMap_ASW_MXL_KINGmat)
    checkIdentical(class(mypcs), "pcair")
    checkEqualsNumeric(mypcs$sum.values, 100.3638, tol=1e-4)
    checkEquals(mypcs$kin.thresh, 2^(-11/2))
    checkEquals(mypcs$div.thresh, -2^(-11/2))
    checkEquals(mypcs$nsamp, 173)
    checkEquals(mypcs$nsnps, 19996)
    checkEquals(mypcs$MAF, 0.01)
    checkEquals(length(mypcs$rels),76)
    checkEquals(length(mypcs$unrels),97)

    obs <- tryCatch( pcair(genoData = HapMap_genoData, kinMat = HapMap_ASW_MXL_KINGmat, divMat = HapMap_ASW_MXL_KINGmat, MAF=-1), error=conditionMessage)
    checkIdentical(obs, "MAF must be in [0,0.5]")

    obs <- tryCatch( pcair(genoData = HapMap_genoData, kinMat = HapMap_ASW_MXL_KINGmat, divMat = HapMap_ASW_MXL_KINGmat, scan.include=1:100), error=conditionMessage)
    checkIdentical(obs, "Not all of the scanID in scan.include are in the provided data")

    newMat <- HapMap_ASW_MXL_KINGmat
    colnames(newMat) <- rownames(newMat) <- NULL
    obs <- tryCatch( pcair(genoData = HapMap_genoData, kinMat = newMat, divMat = newMat), error=conditionMessage)
    checkIdentical(obs, "colnames and rownames of kinMat must be individual IDs")

    obs <- tryCatch( pcair(genoData = HapMap_genoData, kinMat = HapMap_ASW_MXL_KINGmat, divMat = newMat), error=conditionMessage)
    checkIdentical(obs, "colnames and rownames of divMat must be individual IDs")

    #colnames(newMat) <- rownames(newMat) <- 1:ncol(newMat)
    #obs <- tryCatch( pcair(genoData = HapMap_genoData, kinMat = newMat, divMat = newMat), error=conditionMessage)
    #checkIdentical(obs, "")

    obs <- tryCatch( pcair(genoData = HapMap_genoData, kinMat = HapMap_ASW_MXL_KINGmat, divMat = HapMap_ASW_MXL_KINGmat, unrel.set=1:100), error=conditionMessage)
    checkIdentical(obs, "All of the scanIDs in unrel.set must be in the scanIDs of genoData (and scan.include if specified)")

    
    close(HapMap_genoData)
}

