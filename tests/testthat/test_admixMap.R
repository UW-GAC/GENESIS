context("admixMap tests")

test_that("admixMap", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
    gds <- openfn.gds(gdsfile)
    samp <- as.character(read.gdsn(index.gdsn(gds, "sample.id")))
    nsnp <- objdesp.gdsn(index.gdsn(gds, "snp.id"))$dim
    nsamp <- objdesp.gdsn(index.gdsn(gds, "sample.id"))$dim
    closefn.gds(gds)
    set.seed(200); dosage_eur <- sample(0:2, nsnp*nsamp, replace=TRUE)
    set.seed(201); dosage_afr <- ifelse(dosage_eur == 2, 0, sample(0:1, nsnp*nsamp, replace=TRUE))
    set.seed(202); dosage_amer <- 2 - dosage_eur - dosage_afr
    dosage <- list(dosage_eur, dosage_afr, dosage_amer)
    tmpfile <- character(3)
    tmpfile2 <- character(3)
    for (i in 1:3) {
        tmpfile[i] <- tempfile()
        file.copy(gdsfile, tmpfile[i])
        SNPRelate::snpgdsTranspose(tmpfile[i], verbose=FALSE)
        gds <- openfn.gds(tmpfile[i], readonly=FALSE)
        write.gdsn(index.gdsn(gds, "genotype"), matrix(dosage[[i]], nrow=nsamp, ncol=nsnp))
        # factor to character
        delete.gdsn(index.gdsn(gds, "sample.id"))
        add.gdsn(gds, "sample.id", samp)
        add.gdsn(gds, "snp.allele", rep("A,A", nsnp))
        closefn.gds(gds)

        # convert to SeqArray
        tmpfile2[i] <- tempfile()
        seqSNP2GDS(tmpfile[i], tmpfile2[i], verbose=FALSE)
    }

    set.seed(203); pheno <- rnorm(nsamp, mean = 0, sd = 1)
    set.seed(204); covar <- sample(0:1, nsamp, replace=TRUE)
    
    annot <- GWASTools::ScanAnnotationDataFrame(data.frame(scanID = samp, 
                                                    covar, pheno, stringsAsFactors=FALSE))
    genoIterators <- lapply(tmpfile, function(x) {
        gr <- GdsGenotypeReader(x)
        gd <- GenotypeData(gr, scanAnnot=annot)
        GenotypeBlockIterator(gd)
    })
    
    annot <- AnnotatedDataFrame(data.frame(sample.id = samp, 
                                           covar, pheno, stringsAsFactors=FALSE))
    seqIterators <- lapply(tmpfile2, function(x) {
        gr <- seqOpen(x)
        gd <- SeqVarData(gr, sampleData=annot)
        SeqVarBlockIterator(gd, verbose=FALSE)
    })
    
    null.model <- fitNullModel(annot, outcome = "pheno", covars = "covar")
    myassoc <- admixMap(genoIterators, null.model, verbose=FALSE)
    myassoc2 <- admixMap(seqIterators, null.model, verbose=FALSE)
    expect_equal(myassoc, myassoc2)

    # make sure we're reading variant info correctly
    expect_false(any(duplicated(myassoc$variant.id)))

    # check running with only one ancestry
    GWASTools::resetIterator(genoIterators[[1]])
    myassoc3 <- admixMap(genoIterators[[1]], null.model, verbose=FALSE)
    
    lapply(tmpfile, unlink)
    lapply(tmpfile2, unlink)
})
