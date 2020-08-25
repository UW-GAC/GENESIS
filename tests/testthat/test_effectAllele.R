context("effect allele tests")

test_that("effect allele for GenotypeData", {
    genoData <- .testGenoData()
    eff <- effectAllele(genoData)
    expect_true(is.character(eff$effect.allele))
    
    var <- getSnpID(genoData)[1:10]
    eff.var <- effectAllele(genoData, var)
    expect_equal(eff.var, eff[1:10,])

    close(genoData)
})


test_that("effect allele for SeqVarData", {
    svd <- .testData()
    eff <- effectAllele(svd)
    
    var <- seqGetData(svd, "variant.id")[1:10]
    eff.var <- effectAllele(svd, var)
    expect_equal(eff.var, eff[1:10,])
    
    seqClose(svd)
})
