context("Compare with old functions")
library(SeqVarTools)
library(GenomicRanges)
library(Biobase)

## flag for binary tests
## these have unpredictable errors depending on the random covariance matrix
## mostly of the form "a is 0-dimensional"
test_binary <- FALSE

## unpredictable failures - check this
test_GxE <- FALSE

test_that("fitNullModel matches fitNulMM", {
    svd <- .testData()
    grm <- .testGRM(svd)
    gen1 <- GENESIS::fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, verbose=FALSE)
    gen2 <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), cov.mat=grm, verbose=FALSE)

    expect_equivalent(gen2$fixef, gen1$fixef)
    expect_equivalent(gen2$betaCov, gen1$betaCov)
    expect_equivalent(gen2$resid.marginal, gen1$resid.marginal)
    expect_equivalent(gen2$logLik, gen1$logLik)
    expect_equivalent(gen2$logLikR, gen1$logLikR)
    expect_equivalent(gen2$AIC, gen1$AIC)
    expect_equivalent(gen2$workingY, gen1$workingY)
    expect_equivalent(gen2$model.matrix, gen1$model.matrix)
    expect_equivalent(gen2$varComp, gen1$varComp)
    expect_equivalent(gen2$varCompCov, gen1$varCompCov)
    expect_equivalent(gen2$family$family, gen1$family$family)
    expect_equivalent(gen2$zeroFLAG, gen1$zeroFLAG)
    expect_true(all(abs(gen2$cholSigmaInv - gen1$cholSigmaInv) < 1e-9))
    expect_equivalent(gen2$RSS, gen1$RSS)

    seqClose(svd)
})

test_that("fitNullModel matches fitNulMM - group", {
    svd <- .testData()
    grm <- .testGRM(svd)
    gen1 <- GENESIS::fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, group.var="status", verbose=FALSE)
    gen2 <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), cov.mat=grm, group.var="status", verbose=FALSE)

    expect_equivalent(gen2$fixef, gen1$fixef)
    expect_equivalent(gen2$betaCov, gen1$betaCov)
    expect_equivalent(gen2$resid.marginal, gen1$resid.marginal)
    expect_equivalent(gen2$logLik, gen1$logLik)
    expect_equivalent(gen2$logLikR, gen1$logLikR)
    expect_equivalent(gen2$AIC, gen1$AIC)
    expect_equivalent(gen2$workingY, gen1$workingY)
    expect_equivalent(gen2$model.matrix, gen1$model.matrix)
    expect_equivalent(gen2$varComp, gen1$varComp)
    expect_equivalent(gen2$varCompCov, gen1$varCompCov)
    expect_equivalent(gen2$family$family, gen1$family$family)
    expect_equivalent(gen2$zeroFLAG, gen1$zeroFLAG)
    expect_true(all(abs(gen2$cholSigmaInv - gen1$cholSigmaInv) < 1e-9))
    expect_equivalent(gen2$RSS, gen1$RSS)

    seqClose(svd)
})

if (test_binary) {
test_that("fitNullModel matches fitNulMM - binary", {
    svd <- .testData()
    grm <- .testGRM(svd)
    gen1 <- GENESIS::fitNullMM(sampleData(svd), outcome="status", covars=c("sex", "age"), covMatList=grm, family="binomial", verbose=FALSE)
    gen2 <- fitNullModel(svd, outcome="status", covars=c("sex", "age"), cov.mat=grm, family="binomial", verbose=FALSE)

    expect_equivalent(gen2$fixef, gen1$fixef, tolerance=1e-6)
    expect_equivalent(gen2$betaCov, gen1$betaCov, tolerance=1e-6)
    expect_equivalent(gen2$resid.marginal, gen1$resid.marginal, tolerance=1e-6)
    expect_equivalent(gen2$logLik, gen1$logLik, tolerance=1e-6)
    expect_equivalent(gen2$logLikR, gen1$logLikR, tolerance=1e-6)
    expect_equivalent(gen2$AIC, gen1$AIC, tolerance=1e-6)
    expect_equivalent(gen2$workingY, gen1$workingY, tolerance=1e-6)
    expect_equivalent(gen2$model.matrix, gen1$model.matrix, tolerance=1e-6)
    expect_equivalent(gen2$varComp, gen1$varComp, tolerance=1e-6)
    expect_equivalent(gen2$varCompCov, gen1$varCompCov, tolerance=1e-4)
    expect_equivalent(gen2$family$family, gen1$family$family)
    expect_equivalent(gen2$zeroFLAG, gen1$zeroFLAG)
    expect_true(all(abs(gen2$cholSigmaInv - gen1$cholSigmaInv) < 1e-6))
    
    seqClose(svd)
})
}

test_that("fitNullModel matches fitNullReg", {
    svd <- .testData()
    lmm.genesis <- GENESIS::fitNullReg(sampleData(svd), outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)

    expect_true(all(abs(nullmod$fixef - lmm.genesis$fixef) < 1e-9))
    expect_true(all(abs(nullmod$betaCov - lmm.genesis$betaCov) < 1e-9))
    expect_true(all(abs(nullmod$resid.marginal - lmm.genesis$resid.response) < 1e-9))
    expect_true(all(abs(nullmod$logLik - lmm.genesis$logLik) < 1e-9))
    expect_true(all(abs(nullmod$AIC - lmm.genesis$AIC) < 1e-9))
    expect_true(all(nullmod$workingY == lmm.genesis$workingY))
    expect_true(all(nullmod$model.matrix == lmm.genesis$model.matrix))
    expect_equal(nullmod$family$family, lmm.genesis$family$family)

    seqClose(svd)
})

if (test_binary) {
test_that("fitNullModel matches fitNullReg - binary", {
    svd <- .testData()
    lmm.genesis <- GENESIS::fitNullReg(sampleData(svd), outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
    nullmod <- fitNullModel(svd, outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)

    expect_true(all(abs(nullmod$fixef - lmm.genesis$fixef) < 1e-9))
    expect_true(all(abs(nullmod$betaCov - lmm.genesis$betaCov) < 1e-9))
    expect_true(all(abs(nullmod$resid.marginal - lmm.genesis$resid.response) < 1e-9))
    expect_true(all(abs(nullmod$logLik - lmm.genesis$logLik) < 1e-9))
    expect_true(all(abs(nullmod$AIC - lmm.genesis$AIC) < 1e-9))
    expect_true(all(nullmod$workingY == lmm.genesis$workingY))
    expect_true(all(nullmod$model.matrix == lmm.genesis$model.matrix))
    expect_equal(nullmod$family$family, lmm.genesis$family$family)

    seqClose(svd)
})
}


test_that("assocTestSingle matches assocTestMM - Wald", {
    svd <- .testData()
    grm <- .testGRM(svd)
    seqResetFilter(svd, verbose=FALSE)
    nullmod <- GENESIS::fitNullMM(sampleData(svd), outcome="outcome", covMatList=grm, verbose=FALSE)
    
    # multiallelic variants are handled differently
    snv <- isSNV(svd, biallelic=TRUE)
    assoc1 <- GENESIS::assocTestMM(svd, nullmod, test="Wald", snp.include=which(snv), verbose=FALSE)
    
    nullmod <- fitNullModel(svd, outcome="outcome", cov.mat=grm, verbose=FALSE)
    seqSetFilter(svd, variant.sel=snv, verbose=FALSE)
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    assoc2 <- assocTestSingle(iterator, nullmod, test="Wald", verbose=FALSE)
    not1 <- setdiff(assoc2$variant.id, assoc1$snpID)
    expect_true(all(assoc2$freq[not1] %in% c(0,1)))
    expect_true(all(is.na(assoc2$Est[not1])))
    assoc2 <- assoc2[!(assoc2$variant.id %in% not1),]
    expect_equal(nrow(assoc1), nrow(assoc2))
    expect_equal(assoc1$snpID, assoc2$variant.id)
    expect_equal(assoc1$chr, assoc2$chr)
    expect_equal(assoc1$MAF, pmin(assoc2$freq, 1-assoc2$freq))
    expect_equal(assoc1$Est, assoc2$Est)
    expect_equal(assoc1$SE, assoc2$Est.SE)
    expect_equal(assoc1$Wald.Stat, (assoc2$Wald.Stat)^2)
    expect_equal(assoc1$Wald.pval, assoc2$Wald.pval)
    
    seqClose(svd)
})


test_that("assocTestSingle matches assocTestMM - Score", {
    svd <- .testData()
    grm <- .testGRM(svd)
    seqResetFilter(svd, verbose=FALSE)
    nullmod <- GENESIS::fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, verbose=FALSE)
    
    # multiallelic variants are handled differently
    snv <- isSNV(svd, biallelic=TRUE)
    assoc1 <- GENESIS::assocTestMM(svd, nullmod, test="Score", snp.include=which(snv), verbose=FALSE)
    
    nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), cov.mat=grm, verbose=FALSE)
    seqSetFilter(svd, variant.sel=snv, verbose=FALSE)
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    assoc2 <- assocTestSingle(iterator, nullmod, test="Score", verbose=FALSE)
    not1 <- setdiff(assoc2$variant.id, assoc1$snpID)
    assoc2 <- assoc2[!(assoc2$variant.id %in% not1),]
    expect_equal(nrow(assoc1), nrow(assoc2))
    expect_equal(assoc1$Score, assoc2$Score)
    expect_equal(assoc1$Var, (assoc2$Score.SE)^2)
    expect_equal(assoc1$Score.Stat, (assoc2$Score.Stat)^2)
    expect_equal(assoc1$Score.pval, assoc2$Score.pval)
    
    seqClose(svd)
})


if (test_binary) {
test_that("assocTestSingle matches assocTestMM - binary", {
    svd <- .testData()
    grm <- .testGRM(svd)
    seqResetFilter(svd, verbose=FALSE)
    nullmod <- GENESIS::fitNullMM(sampleData(svd), outcome="status", covars=c("sex", "age"), covMatList=grm, family="binomial", verbose=FALSE)
    
    # multiallelic variants are handled differently
    snv <- isSNV(svd, biallelic=TRUE)
    assoc1 <- GENESIS::assocTestMM(svd, nullmod, test="Score", snp.include=which(snv), verbose=FALSE)
    
    nullmod <- fitNullModel(svd, outcome="status", covars=c("sex", "age"), cov.mat=grm, family="binomial", verbose=FALSE)
    seqSetFilter(svd, variant.sel=snv, verbose=FALSE)
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    assoc2 <- assocTestSingle(iterator, nullmod, test="Score", verbose=FALSE)
    not1 <- setdiff(assoc2$variant.id, assoc1$snpID)
    assoc2 <- assoc2[!(assoc2$variant.id %in% not1),]
    expect_equal(nrow(assoc1), nrow(assoc2))
    expect_equal(assoc1$Score, assoc2$Score, tolerance=1e-6)
    expect_equal(assoc1$Var, (assoc2$Score.SE)^2, tolerance=1e-6)
    expect_equal(assoc1$Score.Stat, (assoc2$Score.Stat)^2, tolerance=1e-6)
    expect_equal(assoc1$Score.pval, assoc2$Score.pval, tolerance=1e-6)
    
    seqClose(svd)
})
}


if (test_GxE) {
test_that("assocTestSingle matches assocTestMM - GxE", {
    svd <- .testData()
    grm <- .testGRM(svd)
    seqSetFilter(svd, variant.sel=20:100, verbose=FALSE)
    nullmod <- GENESIS::fitNullMM(sampleData(svd), outcome="outcome", covars="sex", covMatList=grm, verbose=FALSE)
    
    # multiallelic variants are handled differently
    assoc1 <- GENESIS::assocTestMM(svd, nullmod, ivars="sex", verbose=FALSE)
    
    nullmod <- fitNullModel(svd, outcome="outcome", covars="sex", cov.mat=grm, verbose=FALSE)
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    assoc2 <- assocTestSingle(iterator, nullmod, test="Wald", GxE="sex", verbose=FALSE)

    keep <- which(!is.na(assoc1$Est.G) & !is.na(assoc2$Est.G))
    assoc1 <- assoc1[keep,]
    assoc2 <- assoc2[keep,]
    
    expect_equal(assoc1$Est.G, assoc2$Est.G)
    expect_equal(assoc1$`Est.G:sexM`, assoc2$`Est.G:sexM`)
    expect_equal(assoc1$SE.G, assoc2$SE.G)
    expect_equal(assoc1$`SE.G:sexM`, assoc2$`SE.G:sexM`)
    expect_equal(assoc1$GxE.Stat, (assoc2$GxE.Stat)^2)
    expect_equal(assoc1$GxE.pval, assoc2$GxE.pval)
    expect_equal(assoc1$Joint.Stat, (assoc2$Joint.Stat)^2)
    expect_equal(assoc1$Joint.pval, assoc2$Joint.pval)
    
    seqClose(svd)
})
}


test_that("assocTestAggregate matches assocTestSeqWindow - Burden, Wald", {
    svd <- .testData()
    nullmod <- GENESIS::fitNullReg(sampleData(svd), outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc1 <- GENESIS::assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="Burden", burden.test="Wald", verbose=FALSE, chromosome=1)
    
    nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
    assoc2 <- assocTestAggregate(iterator, nullmod, test="Burden", burden.test="Wald", verbose=FALSE)

    var1 <- assoc1$variantInfo[!(assoc1$variantInfo$freq %in% c(0,1)),]
    var2 <- do.call(rbind, assoc2$variantInfo)
    var2 <- var2[!duplicated(paste(var2$variant.id, var2$allele.index)),]
    expect_equal(nrow(var1), nrow(var2))  
    expect_equal(var1$variantID, var2$variant.id)
    expect_equal(var1$allele, var2$allele.index)
    expect_equal(var1$chr, var2$chr)
    expect_equal(var1$pos, var2$pos)
    expect_equal(var1$n.obs, var2$n.obs)
    expect_equal(var1$freq, var2$freq)
    expect_equal(var1$weight, var2$weight)

    res1 <- assoc1$results[assoc1$results$dup %in% 0,]
    res2 <- assoc2$results[assoc2$results$n.site > 0,]
    expect_equal(res1$n.site, res2$n.site)
    expect_equal(res1$Est, res2$Est)
    expect_equal(res1$SE, res2$Est.SE)
    expect_equal(res1$Wald.stat, (res2$Wald.Stat)^2)
    expect_equal(res1$Wald.pval, res2$Wald.pval)
    
    seqClose(svd)
})


test_that("assocTestAggregate matches assocTestSeqWindow - Burden, Score", {
    svd <- .testData()
    nullmod <- GENESIS::fitNullReg(sampleData(svd), outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc1 <- GENESIS::assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="Burden", burden.test="Score", verbose=FALSE, chromosome=1)
    
    nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
    assoc2 <- assocTestAggregate(iterator, nullmod, test="Burden", burden.test="Score", verbose=FALSE)

    res1 <- assoc1$results[assoc1$results$dup %in% 0,]
    res2 <- assoc2$results[assoc2$results$n.site > 0,]
    expect_equal(res1$n.site, res2$n.site)
    expect_equal(res1$Score, res2$Score)
    expect_equal(res1$Var, (res2$Score.SE)^2)
    expect_equal(res1$Score.stat, (res2$Score.Stat)^2)
    expect_equal(res1$Score.pval, res2$Score.pval)
    
    seqClose(svd)
})


test_that("assocTestAggregate matches assocTestSeqWindow - Burden, Wald, LMM", {
    svd <- .testData()
    grm <- .testGRM(svd)
    nullmod <- GENESIS::fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, verbose=FALSE)
    assoc1 <- GENESIS::assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="Burden", burden.test="Wald", verbose=FALSE, chromosome=1)
    
    nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), cov.mat=grm, verbose=FALSE)
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
    assoc2 <- assocTestAggregate(iterator, nullmod, test="Burden", burden.test="Wald", verbose=FALSE)

    res1 <- assoc1$results[assoc1$results$dup %in% 0,]
    res2 <- assoc2$results[assoc2$results$n.site > 0,]
    expect_equal(res1$n.site, res2$n.site)
    expect_equal(res1$Est, res2$Est)
    expect_equal(res1$SE, res2$Est.SE)
    expect_equal(res1$Wald.stat, (res2$Wald.Stat)^2)
    expect_equal(res1$Wald.pval, res2$Wald.pval)
    
    seqClose(svd)
})


test_that("assocTestAggregate matches assocTestSeqWindow - Burden, Score, LMM", {
    svd <- .testData()
    grm <- .testGRM(svd)
    nullmod <- GENESIS::fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, verbose=FALSE)
    assoc1 <- GENESIS::assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="Burden", burden.test="Score", verbose=FALSE, chromosome=1)
    
    nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), cov.mat=grm, verbose=FALSE)
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
    assoc2 <- assocTestAggregate(iterator, nullmod, test="Burden", burden.test="Score", verbose=FALSE)

    res1 <- assoc1$results[assoc1$results$dup %in% 0,]
    res2 <- assoc2$results[assoc2$results$n.site > 0,]
    expect_equal(res1$n.site, res2$n.site)
    expect_equal(res1$Score, res2$Score)
    expect_equal(res1$Var, (res2$Score.SE)^2)
    expect_equal(res1$Score.stat, (res2$Score.Stat)^2)
    expect_equal(res1$Score.pval, res2$Score.pval)
    
    seqClose(svd)
})

test_that("assocTestAggregate matches assocTestSeqWindow - Burden, Binary", {
    svd <- .testData()
    nullmod <- GENESIS::fitNullReg(sampleData(svd), outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
    assoc1 <- GENESIS::assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="Burden", burden.test="Score", verbose=FALSE, chromosome=1)
    
    nullmod <- fitNullModel(svd, outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
    assoc2 <- assocTestAggregate(iterator, nullmod, test="Burden", burden.test="Score", verbose=FALSE)

    res1 <- assoc1$results[assoc1$results$dup %in% 0,]
    res2 <- assoc2$results[assoc2$results$n.site > 0,]
    expect_equal(res1$n.site, res2$n.site)
    expect_equal(res1$Score, res2$Score, tolerance=1e-6)
    expect_equal(res1$Var, (res2$Score.SE)^2, tolerance=1e-6)
    expect_equal(res1$Score.stat, (res2$Score.Stat)^2, tolerance=1e-6)
    expect_equal(res1$Score.pval, res2$Score.pval, tolerance=1e-6)
    
    seqClose(svd)
})


if (test_binary) {
test_that("assocTestAggregate matches assocTestSeqWindow - Burden, Binary, LMM", {
    svd <- .testData()
    grm <- .testGRM(svd)
    nullmod <- GENESIS::fitNullMM(sampleData(svd), outcome="status", covars=c("sex", "age"), covMatList=grm, family="binomial", verbose=FALSE)
    assoc1 <- GENESIS::assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="Burden", burden.test="Score", verbose=FALSE, chromosome=1)
    
    nullmod <- fitNullModel(svd, outcome="status", covars=c("sex", "age"), cov.mat=grm, family="binomial", verbose=FALSE)
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
    assoc2 <- assocTestAggregate(iterator, nullmod, test="Burden", burden.test="Score", verbose=FALSE)

    res1 <- assoc1$results[assoc1$results$dup %in% 0,]
    res2 <- assoc2$results[assoc2$results$n.site > 0,]
    expect_equal(res1$n.site, res2$n.site)
    expect_equal(res1$Score, res2$Score, tolerance=1e-6)
    expect_equal(res1$Var, (res2$Score.SE)^2, tolerance=1e-6)
    expect_equal(res1$Score.stat, (res2$Score.Stat)^2, tolerance=1e-6)
    expect_equal(res1$Score.pval, res2$Score.pval, tolerance=1e-6)
    
    seqClose(svd)
})
}


test_that("assocTestAggregate matches assocTestSeqWindow - SKAT", {
    svd <- .testData()
    nullmod <- GENESIS::fitNullReg(sampleData(svd), outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc1 <- GENESIS::assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="SKAT", verbose=FALSE, chromosome=1)
    
    nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
    assoc2 <- assocTestAggregate(iterator, nullmod, test="SKAT", verbose=FALSE)

    res1 <- assoc1$results[assoc1$results$dup %in% 0,]
    res2 <- assoc2$results[assoc2$results$n.site > 0,]
    expect_equal(res1$n.site, res2$n.site)
    expect_equal(res1$Q_0, res2$Q_0)
    expect_equal(res1$pval_0, res2$pval_0)
    expect_equal(res1$err_0, res2$err_0)
    
    seqClose(svd)
})


if (test_binary) {
test_that("assocTestAggregate matches assocTestSeqWindow - SKAT, binary", {
    svd <- .testData()
    nullmod <- GENESIS::fitNullReg(sampleData(svd), outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
    assoc1 <- GENESIS::assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="SKAT", verbose=FALSE, chromosome=1)
    
    nullmod <- fitNullModel(svd, outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
    assoc2 <- assocTestAggregate(iterator, nullmod, test="SKAT", verbose=FALSE)

    res1 <- assoc1$results[assoc1$results$dup %in% 0,]
    res2 <- assoc2$results[assoc2$results$n.site > 0,]
    expect_equal(res1$n.site, res2$n.site)
    expect_equal(res1$Q_0, res2$Q_0, tolerance=1e-6)
    expect_equal(res1$pval_0, res2$pval_0, tolerance=1e-6)
    expect_equal(res1$err_0, res2$err_0, tolerance=1e-6)
    
    seqClose(svd)
})
}


test_that("assocTestAggregate matches assocTestSeqWindow - SKAT, LMM", {
    svd <- .testData()
    grm <- .testGRM(svd)
    nullmod <- GENESIS::fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, verbose=FALSE)
    assoc1 <- GENESIS::assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="SKAT", verbose=FALSE, chromosome=1)
    
    nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), cov.mat=grm, verbose=FALSE)
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
    assoc2 <- assocTestAggregate(iterator, nullmod, test="SKAT", verbose=FALSE)

    res1 <- assoc1$results[assoc1$results$dup %in% 0,]
    res2 <- assoc2$results[assoc2$results$n.site > 0,]
    expect_equal(res1$n.site, res2$n.site)
    expect_equal(res1$Q_0, res2$Q_0)
    expect_equal(res1$pval_0, res2$pval_0)
    expect_equal(res1$err_0, res2$err_0)
    
    seqClose(svd)
})


if (test_binary) {
test_that("assocTestAggregate matches assocTestSeqWindow - SKAT, binary, LMM", {
    svd <- .testData()
    grm <- .testGRM(svd)
    nullmod <- GENESIS::fitNullMM(sampleData(svd), outcome="status", covars=c("sex", "age"), covMatList=grm, family="binomial", verbose=FALSE)
    assoc1 <- GENESIS::assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="SKAT", verbose=FALSE, chromosome=1)
    
    nullmod <- fitNullModel(svd, outcome="status", covars=c("sex", "age"), cov.mat=grm, family="binomial", verbose=FALSE)
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
    assoc2 <- assocTestAggregate(iterator, nullmod, test="SKAT", verbose=FALSE)

    res1 <- assoc1$results[assoc1$results$dup %in% 0,]
    res2 <- assoc2$results[assoc2$results$n.site > 0,]
    expect_equal(res1$n.site, res2$n.site)
    expect_equal(res1$Q_0, res2$Q_0, tolerance=1e-6)
    expect_equal(res1$pval_0, res2$pval_0, tolerance=1e-6)
    expect_equal(res1$err_0, res2$err_0, tolerance=1e-6)
    
    seqClose(svd)
})
}


test_that("assocTestAggregate matches assocTestSeqWindow - SKAT-O, LMM", {
    svd <- .testData()
    grm <- .testGRM(svd)
    nullmod <- GENESIS::fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, verbose=FALSE)
    assoc1 <- GENESIS::assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="SKAT", rho=c(0,0.5,1), verbose=FALSE, chromosome=1)
    
    nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), cov.mat=grm, verbose=FALSE)
    seqSetFilterChrom(svd, include=1, verbose=FALSE)
    iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
    assoc2 <- assocTestAggregate(iterator, nullmod, test="SKAT", rho=c(0,0.5,1), verbose=FALSE)

    res1 <- assoc1$results[assoc1$results$dup %in% 0,]
    res2 <- assoc2$results[assoc2$results$n.site > 0,]
    expect_equal(res1$n.site, res2$n.site)
    expect_equal(res1$Q_0, res2$Q_0)
    expect_equal(res1$pval_0, res2$pval_0)
    expect_equal(res1$err_0, res2$err_0)
    
    seqClose(svd)
})
