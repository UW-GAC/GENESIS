context("fast score SE")
library(SeqVarTools)

BPPARAM <- BiocParallel::SerialParam()

test_that("fast score SE", {
    svd <- .testData()
    covmat <- .testGRM(svd)
    
    nm <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), cov.mat=covmat, verbose=FALSE)
    set.seed(123)
    nm.se <- fitNullModelFastScore(svd, outcome="outcome", covars=c("sex", "age"), cov.mat=covmat, verbose=FALSE)
    expect_true(all(c("se.correction", "score.table") %in% names(nm.se)))
    
    chk <- intersect(names(nm), names(nm.se))
    expect_equal(nm[chk], nm.se[chk])
    
    set.seed(456)
    score.table <- calcScore(svd, nm, verbose=FALSE)
    nm2 <- nullModelFastScore(nm, score.table, verbose=FALSE)
    expect_equal(nm2, nm.se)
    
    iterator <- SeqVarBlockIterator(svd, verbose=FALSE)
    expect_error(assocTestSingle(iterator, nm, fast.score.SE=TRUE, BPPARAM=BPPARAM, verbose=FALSE),
                 "null.model must have se.correction")
    
    assoc <- assocTestSingle(iterator, nm, BPPARAM=BPPARAM, verbose=FALSE)
    resetIterator(iterator, verbose=FALSE)
    set.seed(789)
    assoc.se <- assocTestSingle(iterator, nm.se, fast.score.SE=TRUE, BPPARAM=BPPARAM, verbose=FALSE)
    expect_equal(assoc, assoc.se, tolerance=0.0001)
    
    seqClose(svd)
})
