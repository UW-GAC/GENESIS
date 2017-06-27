library(GWASdata)
library(SNPRelate)

test_assocTestMM <- function() {
    # use GWASdata objects because they have the Y chromosome)
    data(illuminaScanADF)
    scanAnnot <- illuminaScanADF

    gdsfile <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
    gdsobj <- openfn.gds(gdsfile)
    gds <- GdsGenotypeReader(gdsobj)
    genoData <-  GenotypeData(gds, scanAnnot=scanAnnot)

    grm <- snpgdsGRM(gdsobj, verbose=FALSE)
    covMat <- grm$grm
    dimnames(covMat) <- list(grm$sample.id, grm$sample.id)

    # no Y chrom
    nullmod <- fitNullMM(scanData=scanAnnot, outcome="status", covars="sex", covMatList=covMat, family=binomial, verbose=FALSE)
    assoc <- assocTestMM(genoData, nullmod, test="Score", chromosome=c(21:24,26), verbose=FALSE)

    # Y chrom
    males <- scanAnnot$scanID[scanAnnot$sex == "M"]
    nullmod <- fitNullMM(scanData=scanAnnot, outcome="status", covars=NULL, covMatList=covMat, family=binomial, scan.include=males, verbose=FALSE)
    assoc <- assocTestMM(genoData, nullmod, test="Score", chromosome=25, verbose=FALSE)
    
    close(genoData)
}
