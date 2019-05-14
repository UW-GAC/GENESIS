
.testGenoData <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    
    # use GWASdata objects because they have X and Y chromosomes
    scanAnnot <- get(data("illuminaScanADF", package="GWASdata", envir=environment()))
    set.seed(45); scanAnnot$outcome <- rnorm(nrow(scanAnnot))

    gdsfile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
    gds <- GWASTools::GdsGenotypeReader(gdsfile)
    GWASTools::GenotypeData(gds, scanAnnot=scanAnnot)
}


.testGenoDataGRM <- function(genoData) {
    gdsobj <- genoData@data@handler
    grm <- suppressMessages(SNPRelate::snpgdsGRM(gdsobj, verbose=FALSE))
    covMat <- grm$grm
    dimnames(covMat) <- list(grm$sample.id, grm$sample.id)
    covMat
}


.testGenoPCs <- function(genoData) {
    kinship <- .testKing(genoData)
    mypcair <- suppressWarnings(pcair(genoData, kinobj=kinship, divobj=kinship, verbose=FALSE))
    mypcair$vectors[,1:2]
}


.testHMData <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
    HapMap_geno <- GWASTools::GdsGenotypeReader(gdsfile)
    GWASTools::GenotypeData(HapMap_geno)
}


.testHMPCs <- function(genoData) {
    KINGmat <- get(data("HapMap_ASW_MXL_KINGmat", package="GENESIS", envir=environment()))
    mypcair <- pcair(genoData, kinobj = KINGmat, divobj = KINGmat, verbose=FALSE)
    mypcair$vectors[,1:2]
}
