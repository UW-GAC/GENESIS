
.testData <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    sample.id <- SeqArray::seqGetData(gds, "sample.id")
    df <- data.frame(sample.id=sample.id,
                     sex=sample(c("M","F"), replace=TRUE, length(sample.id)),
                     age=rnorm(length(sample.id), mean=50, sd=10),
                     outcome=rnorm(length(sample.id), mean=10, sd=5),
                     status=rbinom(length(sample.id), size=1, prob=0.4),
                     stringsAsFactors=FALSE)
    SeqVarTools::sampleData(gds) <- Biobase::AnnotatedDataFrame(df)
    gds
}

.testKing <- function(gds){
    if (is(gds, "GenotypeData")) {
        gds <- gds@data@handler
    }
    suppressMessages(ibd <- SNPRelate::snpgdsIBDKING(gds, verbose=FALSE))
    kc <- ibd$kinship
    rownames(kc) <- ibd$sample.id
    colnames(kc) <- ibd$sample.id
    kc
}

.testGRM <- function(seqData, ...){
    kinship <- .testKing(seqData)
    mypcair <- suppressWarnings(pcair(seqData, kinobj=kinship, divobj=kinship, verbose=FALSE, ...))
    mypcrel <- pcrelate(seqData, pcMat=mypcair$vectors[,1:2], training.set=mypcair$unrels, verbose=FALSE, ...)
    pcrelateMakeGRM(mypcrel)
}

.testPCs <- function(seqData, ...){
    kinship <- .testKing(seqData)
    mypcair <- suppressWarnings(pcair(seqData, kinobj=kinship, divobj=kinship, verbose=FALSE, ...))
    mypcair$vectors[,1:2]
}
