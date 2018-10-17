context("SeqVarData tests")
library(SeqVarTools)
library(GWASTools)
library(Biobase)
library(gdsfmt)

.testObjects <- function(){
    showfile.gds(closeall=TRUE, verbose=FALSE)
    
    gdsfile <- seqExampleFileName("gds")
    gds.seq <- seqOpen(gdsfile)
    seqSetFilter(gds.seq, variant.sel=isSNV(gds.seq), verbose=FALSE)
    
    data(pedigree)
    pedigree <- pedigree[match(seqGetData(gds.seq, "sample.id"), pedigree$sample.id),]
    pedigree$outcome <- rnorm(nrow(pedigree))
    seqData <- SeqVarData(gds.seq, sampleData=AnnotatedDataFrame(pedigree))
    
    snpfile <- tempfile()
    seqGDS2SNP(gds.seq, snpfile, verbose=FALSE)
    
    gds.snp <- GdsGenotypeReader(snpfile)  
    names(pedigree)[names(pedigree) == "sample.id"] <- "scanID"
    genoData <- GenotypeData(gds.snp, scanAnnot=ScanAnnotationDataFrame(pedigree))

    list(seqData=seqData, genoData=genoData)
}

.testCleanup <- function(data){
    seqClose(data$seqData)
    close(data$genoData)
    unlink(data$genoData@data@filename)
}

.testPCair <- function(genoData, kinship, ...){
    pcair(genoData, kinMat=kinship, divMat=kinship, verbose=FALSE, ...)
}

.testPCrelate <- function(genoData, mypcair, ...){
    mypcrel <- pcrelate(genoData, pcMat=mypcair$vectors[,1:2], training.set=mypcair$unrels, verbose=FALSE, ...)
    pcrelateMakeGRM(mypcrel)
}

.testNullModel <- function(genoData, grm){
    if (is(genoData, "GenotypeData")) {
        scanData <- getScanAnnotation(genoData)
    } else if (is(genoData, "SeqVarData")) {
        scanData <- sampleData(genoData)
    }
    fitNullMM(scanData, outcome="outcome", covars="sex", covMatList=grm, verbose=FALSE)
}

test_that("PCair", {
    data <- .testObjects()
    kinship <- .testKing(data$genoData)
    
    pca.seq <- .testPCair(data$seqData, kinship)
    pca.snp <- .testPCair(data$genoData, kinship)
    expect_equal(pca.seq$values, pca.snp$values)
    expect_equal(pca.seq$unrels, pca.snp$unrels)

    snpID <- sample(getSnpID(data$genoData), 100)
    scanID <- sample(getScanID(data$genoData), 50)
    seqSetFilter(data$seqData, verbose=FALSE)
    pca.seq <- .testPCair(data$seqData, kinship, scan.include=scanID, snp.include=snpID)
    pca.snp <- .testPCair(data$genoData, kinship, scan.include=scanID, snp.include=snpID)
    expect_equal(pca.seq$values, pca.snp$values)
    expect_equal(pca.seq$unrels, pca.snp$unrels)
    
    .testCleanup(data)
})

test_that("PCrelate", {
    data <- .testObjects()
    kinship <- .testKing(data$genoData)
    mypcair <- .testPCair(data$genoData, kinship)
    
    grm.seq <- .testPCrelate(data$seqData, mypcair)
    grm.snp <- .testPCrelate(data$genoData, mypcair)
    names(dimnames(grm.seq)) <- NULL
    expect_equal(grm.seq, grm.snp)

    snpID <- sample(getSnpID(data$genoData), 100)
    scanID <- sample(getScanID(data$genoData), 50)
    mypcair <- .testPCair(data$genoData, kinship, scan.include=scanID, snp.include=snpID)
    seqSetFilter(data$seqData, verbose=FALSE)
    grm.seq <- .testPCrelate(data$seqData, mypcair, scan.include=scanID, snp.include=snpID)
    grm.snp <- .testPCrelate(data$genoData, mypcair, scan.include=scanID, snp.include=snpID)
    names(dimnames(grm.seq)) <- NULL
    expect_equal(grm.seq, grm.snp)
    
    .testCleanup(data)
})

test_that("nullModel", {
    data <- .testObjects()
    kinship <- .testKing(data$genoData)
    mypcair <- .testPCair(data$genoData, kinship)
    grm <- .testPCrelate(data$genoData, mypcair)

    null.seq <- .testNullModel(data$seqData, grm)
    null.snp <- .testNullModel(data$genoData, grm)
    expect_equal(null.seq$varComp, null.snp$varComp)
    
    .testCleanup(data)
})

test_that("getSnpIndex", {
    data <- .testObjects()

    snp.seq <- getSnpIndex(data$seqData, snp.include=1:10, chromosome=NULL)
    expect_equal(1:10, snp.seq$index)

    seqSetFilter(data$seqData, variant.id=1:10, verbose=FALSE)
    snp.seq <- getSnpIndex(data$seqData, snp.include=NULL, chromosome=NULL)
    expect_equal(1:10, snp.seq$value)
    expect_equal(1:10, snp.seq$index)
    
    seqSetFilter(data$seqData, variant.id=11:20, verbose=FALSE)
    snp.seq <- getSnpIndex(data$seqData, snp.include=NULL, chromosome=NULL)
    expect_equal(11:20, snp.seq$value)
    expect_equal(1:10, snp.seq$index)
    
    seqSetFilter(data$seqData, verbose=FALSE)
    snp.seq <- getSnpIndex(data$seqData, snp.include=NULL, chromosome=22)
    chr22 <- which(seqGetData(data$seqData, "chromosome") == 22)
    expect_equal(chr22, snp.seq$index)    
    
    .testCleanup(data)
})

test_that("prepareGenotype", {
    data <- .testObjects()

    snp.id <- sample(getSnpID(data$genoData), 50)
    scan.id <- sample(getScanID(data$genoData), 20)
    seqResetFilter(data$seqData, verbose=FALSE)
    snp.idx <- which(seqGetData(data$seqData, "variant.id") %in% snp.id)
    scan.idx <- which(seqGetData(data$seqData, "sample.id") %in% scan.id)
    geno.seq <- prepareGenotype(data$seqData, snp.idx, scan.idx, impute.geno=TRUE)
    snp.idx <- which(getSnpID(data$genoData) %in% snp.id)
    scan.idx <- which(getScanID(data$genoData) %in% scan.id)
    geno.snp <- prepareGenotype(data$genoData, snp.idx, scan.idx, impute.geno=TRUE)
    names(dimnames(geno.seq)) <- NULL
    expect_equal(geno.seq, 2-geno.snp)

    .testCleanup(data)
})

test_that("assocMM", {
    data <- .testObjects()

    kinship <- .testKing(data$genoData)
    mypcair <- .testPCair(data$genoData, kinship)
    grm <- .testPCrelate(data$genoData, mypcair)
    nullmod <- .testNullModel(data$genoData, grm)

    assoc.seq <- assocTestMM(data$seqData, nullMMobj=nullmod, verbose=FALSE)
    assoc.snp <- assocTestMM(data$genoData, nullMMobj=nullmod, verbose=FALSE)
    cols <- !(names(assoc.snp) %in% c("minor.allele", "Est"))
    expect_equal(assoc.seq[,cols], assoc.snp[,cols])
    expect_equal(assoc.seq$Est, -assoc.snp$Est)
    ind <- assoc.seq$MAF < 0.5
    expect_equal(assoc.seq$minor.allele[ind] == "ref", assoc.snp$minor.allele[ind] == "A")

    snpID <- sample(getSnpID(data$genoData), 100)
    assoc.seq <- assocTestMM(data$seqData, nullMMobj=nullmod, snp.include=snpID, verbose=FALSE)
    assoc.snp <- assocTestMM(data$genoData, nullMMobj=nullmod, snp.include=snpID, verbose=FALSE)
    cols <- !(names(assoc.snp) %in% c("minor.allele", "Est"))
    expect_equal(assoc.seq[,cols], assoc.snp[,cols])
    expect_equal(assoc.seq$Est, -assoc.snp$Est)

    # check filter
    seqSetFilter(data$seqData, variant.id=snpID, verbose=FALSE)
    assoc.seq <- assocTestMM(data$seqData, nullMMobj=nullmod, verbose=FALSE)
    expect_equal(assoc.seq[,cols], assoc.snp[,cols])
    expect_equal(assoc.seq$Est, -assoc.snp$Est)

    .testCleanup(data)
})
