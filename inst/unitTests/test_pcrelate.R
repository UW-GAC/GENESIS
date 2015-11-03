test_pcrelate <- function(){

	# file path to GDS file
    gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
    # read in GDS data
    HapMap_geno <- GdsGenotypeReader(filename = gdsfile)
    # create a GenotypeData class object
    HapMap_genoData <- GenotypeData(HapMap_geno)
    # load saved matrix of KING-robust estimates
    data("HapMap_ASW_MXL_KINGmat")
    # PC-AiR
    mypcs <- pcair(genoData = HapMap_genoData, kinMat = HapMap_ASW_MXL_KINGmat, divMat = HapMap_ASW_MXL_KINGmat)

    myrel <- pcrelate(genoData = HapMap_genoData, pcMat = mypcs$vectors[,1:2])

    close(HapMap_genoData)	
}