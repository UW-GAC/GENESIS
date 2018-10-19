
.testGenoData <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    
    # use GWASdata objects because they have X and Y chromosomes
    scanAnnot <- get(data("illuminaScanADF", package="GWASdata", envir=environment()))
    scanAnnot$outcome <- rnorm(nrow(scanAnnot))

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


.testHMData <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
    HapMap_geno <- GWASTools::GdsGenotypeReader(gdsfile)
    GWASTools::GenotypeData(HapMap_geno)
}


.testHMPCs <- function(genoData) {
    KINGmat <- get(data("HapMap_ASW_MXL_KINGmat", package="GENESIS", envir=environment()))
    mypcs <- pcair(genoData, kinobj = KINGmat, divobj = KINGmat, verbose=FALSE)
    mypcs$vectors[,1:2]
}
