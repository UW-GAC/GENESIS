context("use imputed dosage")
library(SeqVarTools)

.testImputedData <- function() {
    gds <- SeqVarTools:::.testDosageData()
    gds <- SeqVarData(gds)
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

test_that("assocTestSingle - imputed", {
    svd <- .testImputedData()
    seqSetFilter(svd, variant.sel=isSNV(svd), verbose=FALSE)
    iterator <- SeqVarBlockIterator(svd, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
    resetIterator(iterator, verbose=FALSE)
    assoc2 <- assocTestSingle(iterator, nullmod, imputed=TRUE, verbose=FALSE)
    expect_equal(assoc$variant.id, assoc2$variant.id)
    seqClose(svd)
})

test_that("assocTestAggregate - imputed", {
    svd <- .testImputedData()
    seqSetFilter(svd, variant.sel=isSNV(svd), verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=1e5, verbose=FALSE)
    nullmod <- fitNullModel(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestAggregate(iterator, nullmod, verbose=FALSE)
    resetIterator(iterator, verbose=FALSE)
    assoc2 <- assocTestAggregate(iterator, nullmod, imputed=TRUE, verbose=FALSE)
    expect_equal(assoc$variantInfo$variant.id, assoc2$variantInfo$variant.id)
    seqClose(svd)
})
