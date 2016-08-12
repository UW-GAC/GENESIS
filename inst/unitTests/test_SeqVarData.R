library(SeqVarTools)
library(Biobase)
library(SNPRelate)

.testObjects <- function(){
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    
    gdsfile <- seqExampleFileName("gds")
    gds.seq <- seqOpen(gdsfile)
    seqSetFilter(gds.seq, variant.sel=isSNV(gds.seq), verbose=FALSE)
    
    snpfile <- tempfile()
    seqGDS2SNP(gds.seq, snpfile, verbose=FALSE)
    
    data(pedigree)
    gds.seq <- seqOpen(gdsfile)
    pedigree <- pedigree[match(seqGetData(gds.seq, "sample.id"), pedigree$sample.id),]
    pedigree$outcome <- rnorm(nrow(pedigree))
    seqSetFilter(gds.seq, variant.sel=isSNV(gds.seq), verbose=FALSE)
    seqData <- SeqVarData(gds.seq, sampleData=AnnotatedDataFrame(pedigree))
    
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

.testPCair <- function(genoData, kinship){
    pcair(genoData=genoData, kinMat=kinship, divMat=kinship, verbose=FALSE)
}

.testPCrelate <- function(genoData, mypcair){
    mypcrel <- pcrelate(genoData, pcMat=mypcair$vectors[,1:2], training.set=mypcair$unrels, verbose=FALSE)
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

test_PCrelate <- function(){
    data <- .testObjects()
    kinship <- .testKing(data$genoData)
    mypcair <- .testPCair(data$genoData, kinship)
    
    grm.seq <- .testPCrelate(data$seqData, mypcair)
    grm.snp <- .testPCrelate(data$genoData, mypcair)
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

    snp.id <- c(10:20, 30:50)
    scan.id <- getScanID(data$genoData)[c(20:30, 40:50)]
    geno.seq <- GENESIS:::prepareGenotype(data$seqData, snp.id, scan.id, impute.geno=TRUE)
    geno.snp <- GENESIS:::prepareGenotype(data$genoData, snp.id, scan.id, impute.geno=TRUE)
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

    .testCleanup(data)
}
