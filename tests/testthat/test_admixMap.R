context("admixMap tests")
library(SeqVarTools)

.testLocal <- function(type=c("GenotypeData", "SeqVarData")) {
    type <- match.arg(type)
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
    gds <- openfn.gds(gdsfile)
    samp <- read.gdsn(index.gdsn(gds, "sample.id"))
    nsnp <- objdesp.gdsn(index.gdsn(gds, "snp.id"))$dim
    nsamp <- objdesp.gdsn(index.gdsn(gds, "sample.id"))$dim
    closefn.gds(gds)
    dosage_eur <- sample(0:2, nsnp*nsamp, replace=TRUE)
    dosage_afr <- ifelse(dosage_eur == 2, 0, sample(0:1, nsnp*nsamp, replace=TRUE))
    dosage_amer <- 2 - dosage_eur - dosage_afr
    dosage <- list(dosage_eur, dosage_afr, dosage_amer)
    tmpfile <- character(3)
    for (i in 1:3) {
        tmpfile[i] <- tempfile()
        file.copy(gdsfile, tmpfile[i])
        gds <- openfn.gds(tmpfile[i], readonly=FALSE)
        write.gdsn(index.gdsn(gds, "genotype"), matrix(dosage[[i]], nrow=nsamp, ncol=nsnp))
        closefn.gds(gds)
        if (type == "SeqVarData") {
            tmpfile2 <- tempfile()
            seqSNP2GDS(tmpfile[i], tmpfile2)
            unlink(tmpfile[i])
            tmpfile[i] <- tmpfile2
        }
    }

    pheno <- rnorm(nsamp, mean = 0, sd = 1)
    covar <- sample(0:1, nsamp, replace=TRUE)
    
    if (type == "GenotypeData") {
        annot <- ScanAnnotationDataFrame(data.frame(scanID = samp, 
                                                    covar, pheno, stringsAsFactors=FALSE))
        genoIterators <- lapply(tmpfile, function(x) {
            gr <- GdsGenotypeReader(x)
            gd <- GenotypeData(gr, scanAnnot=annot)
            GenotypeBlockIterator(gd)
        })
    } else {
        annot <- AnnotatedDataFrame(data.frame(sample.id = samp, 
                                               covar, pheno, stringsAsFactors=FALSE))
        genoIterators <- lapply(tmpfile, function(x) {
            gr <- seqOpen(x)
            gd <- GenotypeData(gr, scanAnnot=annot)
            GenotypeBlockIterator(gd)
        })
    }
    setNames(genoIterators, tmpfile)
}


.closeLocal <- function(tmpfile) {
    lapply(tmpfile, unlink)
}


test_that("admixMap - GenotypeData", {
    genoIterators <- .testLocal("GenotypeData")
    null.model <- fitNullModel(genoIterators[[1]], outcome = "pheno", covars = "covar")
    myassoc <- admixMap(genoIterators, null.model)
})
