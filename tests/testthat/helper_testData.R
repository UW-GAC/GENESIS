
.testData <- function() {
    gds <- SeqVarTools:::.testSeqVarData()
    sample.id <- SeqArray::seqGetData(gds, "sample.id")
    set.seed(1); sex <- sample(c("M","F"), replace=TRUE, length(sample.id))
    set.seed(2); age <- rnorm(length(sample.id), mean=50, sd=10)
    set.seed(3); outcome <- rnorm(length(sample.id), mean=10, sd=5)
    set.seed(4); status <- rbinom(length(sample.id), size=1, prob=0.4)
    df <- data.frame(sample.id, sex, age, outcome, status,
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
    seqData <- SeqVarTools:::SeqVarBlockIterator(seqData, verbose=FALSE)
    mypcrel <- pcrelate(seqData, pcs=mypcair$vectors[,1:2], training.set=mypcair$unrels, verbose=FALSE, ...)
    pcrelateToMatrix(mypcrel, verbose=FALSE)
}

.testPCs <- function(seqData, ...){
    kinship <- .testKing(seqData)
    mypcair <- suppressWarnings(pcair(seqData, kinobj=kinship, divobj=kinship, verbose=FALSE, ...))
    mypcair$vectors[,1:2]
}
