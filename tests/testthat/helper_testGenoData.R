
.testGenoData <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    
    # use GWASdata objects because they have X and Y chromosomes
    scanAnnot <- get(data("illuminaScanADF", package="GWASdata", envir=environment()))
    scanAnnot$outcome <- rnorm(nrow(scanAnnot))

    gdsfile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
    gds <- GdsGenotypeReader(gdsfile)
    GenotypeData(gds, scanAnnot=scanAnnot)
}


.testGenoDataGRM <- function(genoData) {
    gdsobj <- genoData@data@handler
    grm <- suppressMessages(SNPRelate::snpgdsGRM(gdsobj, verbose=FALSE))
    covMat <- grm$grm
    dimnames(covMat) <- list(grm$sample.id, grm$sample.id)
    covMat
}
