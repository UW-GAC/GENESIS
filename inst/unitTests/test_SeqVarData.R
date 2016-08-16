library(SeqVarTools)
library(GWASTools)
library(Biobase)
library(SNPRelate)
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

.testKing <- function(gds){
    if (is(gds, "GenotypeData")) {
        gds <- gds@data@handler
    }
    suppressMessages(ibd <- snpgdsIBDKING(gds, verbose=FALSE))
    kc <- ibd$kinship
    rownames(kc) <- ibd$sample.id
    colnames(kc) <- ibd$sample.id
    kc
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

test_PCair <- function(){
    data <- .testObjects()
    kinship <- .testKing(data$genoData)
    
    pca.seq <- .testPCair(data$seqData, kinship)
    pca.snp <- .testPCair(data$genoData, kinship)
    checkEquals(pca.seq$values, pca.snp$values)
    checkEquals(pca.seq$unrels, pca.snp$unrels)

    snpID <- sample(getSnpID(data$genoData), 100)
    scanID <- sample(getScanID(data$genoData), 50)
    seqSetFilter(data$seqData, verbose=FALSE)
    pca.seq <- .testPCair(data$seqData, kinship, scan.include=scanID, snp.include=snpID)
    pca.snp <- .testPCair(data$genoData, kinship, scan.include=scanID, snp.include=snpID)
    checkEquals(pca.seq$values, pca.snp$values)
    checkEquals(pca.seq$unrels, pca.snp$unrels)
    
    .testCleanup(data)
}

test_PCrelate <- function(){
    data <- .testObjects()
    kinship <- .testKing(data$genoData)
    mypcair <- .testPCair(data$genoData, kinship)
    
    grm.seq <- .testPCrelate(data$seqData, mypcair)
    grm.snp <- .testPCrelate(data$genoData, mypcair)
    names(dimnames(grm.seq)) <- NULL
    checkEquals(grm.seq, grm.snp)

    snpID <- sample(getSnpID(data$genoData), 100)
    scanID <- sample(getScanID(data$genoData), 50)
    mypcair <- .testPCair(data$genoData, kinship, scan.include=scanID, snp.include=snpID)
    seqSetFilter(data$seqData, verbose=FALSE)
    grm.seq <- .testPCrelate(data$seqData, mypcair, scan.include=scanID, snp.include=snpID)
    grm.snp <- .testPCrelate(data$genoData, mypcair, scan.include=scanID, snp.include=snpID)
    names(dimnames(grm.seq)) <- NULL
    checkEquals(grm.seq, grm.snp)
    
    .testCleanup(data)
}

test_nullModel <- function(){
    data <- .testObjects()
    kinship <- .testKing(data$genoData)
    mypcair <- .testPCair(data$genoData, kinship)
    grm <- .testPCrelate(data$genoData, mypcair)

    null.seq <- .testNullModel(data$seqData, grm)
    null.snp <- .testNullModel(data$genoData, grm)
    checkEquals(null.seq$varComp, null.snp$varComp)
    
    .testCleanup(data)
}

test_getSnpIndex <- function(){
    data <- .testObjects()

    snp.seq <- GENESIS:::getSnpIndex(data$seqData, snp.include=1:10, chromosome=NULL)
    checkEquals(1:10, snp.seq$index)

    seqSetFilter(data$seqData, variant.id=1:10, verbose=FALSE)
    snp.seq <- GENESIS:::getSnpIndex(data$seqData, snp.include=NULL, chromosome=NULL)
    checkEquals(1:10, snp.seq$value)
    checkEquals(1:10, snp.seq$index)
    
    seqSetFilter(data$seqData, variant.id=11:20, verbose=FALSE)
    snp.seq <- GENESIS:::getSnpIndex(data$seqData, snp.include=NULL, chromosome=NULL)
    checkEquals(11:20, snp.seq$value)
    checkEquals(1:10, snp.seq$index)
    
    seqSetFilter(data$seqData, verbose=FALSE)
    snp.seq <- GENESIS:::getSnpIndex(data$seqData, snp.include=NULL, chromosome=22)
    chr22 <- which(seqGetData(data$seqData, "chromosome") == 22)
    checkEquals(chr22, snp.seq$index)    
    
    .testCleanup(data)
}

test_prepareGenotype <- function(){
    data <- .testObjects()

    snp.id <- sample(getSnpID(data$genoData), 50)
    scan.id <- sample(getScanID(data$genoData), 20)
    seqResetFilter(data$seqData)
    snp.idx <- which(seqGetData(data$seqData, "variant.id") %in% snp.id)
    scan.idx <- which(seqGetData(data$seqData, "sample.id") %in% scan.id)
    geno.seq <- GENESIS:::prepareGenotype(data$seqData, snp.idx, scan.idx, impute.geno=TRUE)
    snp.idx <- which(getSnpID(data$genoData) %in% snp.id)
    scan.idx <- which(getScanID(data$genoData) %in% scan.id)
    geno.snp <- GENESIS:::prepareGenotype(data$genoData, snp.idx, scan.idx, impute.geno=TRUE)
    names(dimnames(geno.seq)) <- NULL
    checkEquals(geno.seq, 2-geno.snp)

    .testCleanup(data)
}

test_assocMM <- function(){
    data <- .testObjects()

    kinship <- .testKing(data$genoData)
    mypcair <- .testPCair(data$genoData, kinship)
    grm <- .testPCrelate(data$genoData, mypcair)
    nullmod <- .testNullModel(data$genoData, grm)

    assoc.seq <- assocTestMM(data$seqData, nullMMobj=nullmod, verbose=FALSE)
    assoc.snp <- assocTestMM(data$genoData, nullMMobj=nullmod, verbose=FALSE)
    cols <- !(names(assoc.snp) %in% c("minor.allele", "Est"))
    checkEquals(assoc.seq[,cols], assoc.snp[,cols])
    checkEquals(assoc.seq$Est, -assoc.snp$Est)
    ind <- assoc.seq$MAF < 0.5
    checkEquals(assoc.seq$minor.allele[ind] == "ref", assoc.snp$minor.allele[ind] == "A")

    snpID <- sample(getSnpID(data$genoData), 100)
    assoc.seq <- assocTestMM(data$seqData, nullMMobj=nullmod, snp.include=snpID, verbose=FALSE)
    assoc.snp <- assocTestMM(data$genoData, nullMMobj=nullmod, snp.include=snpID, verbose=FALSE)
    cols <- !(names(assoc.snp) %in% c("minor.allele", "Est"))
    checkEquals(assoc.seq[,cols], assoc.snp[,cols])
    checkEquals(assoc.seq$Est, -assoc.snp$Est)

    # check filter
    seqSetFilter(data$seqData, variant.id=snpID)
    assoc.seq <- assocTestMM(data$seqData, nullMMobj=nullmod, verbose=FALSE)
    checkEquals(assoc.seq[,cols], assoc.snp[,cols])
    checkEquals(assoc.seq$Est, -assoc.snp$Est)

    .testCleanup(data)
}
