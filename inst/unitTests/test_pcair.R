test_pcair <- function(){
    # file path to GDS file
    gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
    # read in GDS data
    HapMap_geno <- GdsGenotypeReader(filename = gdsfile)
    # create a GenotypeData class object
    HapMap_genoData <- GenotypeData(HapMap_geno)
    # load saved matrix of KING-robust estimates
    data("HapMap_ASW_MXL_KINGmat")
    
    newMat <- HapMap_ASW_MXL_KINGmat
    colnames(newMat) <- NULL
    rownames(newMat) <- NULL
    
    # run PC-AiR
    obs <- tryCatch( pcair(genoData = HapMap_genoData, kinMat = newMat, divMat = newMat), error=conditionMessage)
    checkIdentical("colnames and rownames of kinMat must be individual IDs", obs)
    
    close(HapMap_genoData)
}