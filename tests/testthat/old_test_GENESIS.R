context("Compare with old functions")
library(SeqVarTools)
library(GenomicRanges)
library(Biobase)
library(gdsfmt)

## flags for running various tests

## these have unpredictable errors depending on the random covariance matrix
## mostly of the form "a is 0-dimensional"
test_binary <- TRUE

## unpredictable failures - check this
test_GxE <- TRUE

test_that("fitNullMod matches fitNullReg - linear", {
    dat <- .testNullInputs()
    nullmod <- .fitNullModel(dat$y, dat$X, verbose=FALSE)

    df <- .asDataFrame(dat)
    lm.genesis <- fitNullReg(df, outcome="y", covars=c("X1", "X2", "X3"), verbose=FALSE)

    expect_equivalent(nullmod$fixef, lm.genesis$fixef)
    expect_equivalent(nullmod$betaCov, lm.genesis$betaCov)
    expect_equivalent(nullmod$resid.marginal, lm.genesis$resid.response)
    expect_equivalent(nullmod$logLik, lm.genesis$logLik)
    expect_equivalent(nullmod$AIC, lm.genesis$AIC)
    expect_equivalent(nullmod$workingY, lm.genesis$workingY)
    expect_equivalent(nullmod$model.matrix, lm.genesis$model.matrix)
    expect_equivalent(nullmod$varComp, lm.genesis$sigma^2)
    expect_equal(nullmod$family$family, lm.genesis$family$family)
})


if (test_binary) {
test_that("fitNullMod matches fitNullReg - binary", {
    dat <- .testNullInputs(binary=TRUE)
    nullmod <- .fitNullModel(dat$y, dat$X, family="binomial", verbose=FALSE)

    df <- .asDataFrame(dat)
    glm.genesis <- fitNullReg(df, outcome="y", covars=c("X1", "X2", "X3"), family="binomial", verbose=FALSE)
    
    expect_equivalent(nullmod$fixef, glm.genesis$fixef)
    expect_equivalent(nullmod$betaCov, glm.genesis$betaCov)
    expect_equivalent(nullmod$resid.marginal, glm.genesis$resid.response)
    expect_equivalent(nullmod$logLik, glm.genesis$logLik)
    expect_equivalent(nullmod$AIC, glm.genesis$AIC)
    expect_equivalent(nullmod$workingY, glm.genesis$workingY)
    expect_equivalent(nullmod$model.matrix, glm.genesis$model.matrix)
    expect_equal(nullmod$family$family, glm.genesis$family$family)
    #expect_equivalent(nullmod$varComp, glm.genesis$sigma^2)
    expect_equivalent(diag(as.matrix(nullmod$cholSigmaInv)), glm.genesis$sigma)
})
}


test_that("fitNullMod matches fitNullMM - linear, no group", {
    dat <- .testNullInputs()
    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)

    df <- .asDataFrame(dat)
    lmm.genesis <- fitNullMM(df, outcome="y", covars=c("X1", "X2", "X3"), covMatList=dat$cor.mat, verbose=FALSE)

    expect_equivalent(nullmod$fixef, lmm.genesis$fixef)
    expect_equivalent(nullmod$betaCov, lmm.genesis$betaCov)
    expect_equivalent(nullmod$resid.marginal, lmm.genesis$resid.marginal)
    expect_equivalent(nullmod$logLik, lmm.genesis$logLik)
    expect_equivalent(nullmod$logLikR, lmm.genesis$logLikR)

    expect_equivalent(nullmod$AIC, lmm.genesis$AIC)
    expect_equivalent(nullmod$workingY, lmm.genesis$workingY)
    expect_equivalent(nullmod$model.matrix, lmm.genesis$model.matrix)
    expect_equivalent(nullmod$varComp, lmm.genesis$varComp)
    expect_equivalent(nullmod$varCompCov, lmm.genesis$varCompCov)
    expect_equivalent(nullmod$family$family, lmm.genesis$family$family)
    expect_equivalent(nullmod$zeroFLAG, lmm.genesis$zeroFLAG)
    expect_true(all(abs(nullmod$cholSigmaInv - lmm.genesis$cholSigmaInv) < 1e-7))
    expect_equivalent(nullmod$RSS, lmm.genesis$RSS)
})


test_that("fitNullMod matches fitNullMM - linear, with group", {
    dat <- .testNullInputs()
    nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, group.idx=dat$group.idx, verbose=FALSE)

    df <- .asDataFrame(dat)
    lmm.genesis <- fitNullMM(df, outcome="y", covars=c("X1", "X2", "X3"), covMatList=dat$cor.mat, group.var="group", verbose=FALSE)

    expect_equivalent(nullmod$fixef, lmm.genesis$fixef)
    expect_equivalent(nullmod$betaCov, lmm.genesis$betaCov)
    expect_equivalent(nullmod$resid.marginal, lmm.genesis$resid.marginal)
    expect_equivalent(nullmod$logLik, lmm.genesis$logLik)
    expect_equivalent(nullmod$logLikR, lmm.genesis$logLikR)

    expect_equivalent(nullmod$AIC, lmm.genesis$AIC)
    expect_equivalent(nullmod$workingY, lmm.genesis$workingY)
    expect_equivalent(nullmod$model.matrix, lmm.genesis$model.matrix)
    expect_equivalent(nullmod$varComp, lmm.genesis$varComp)
    expect_equivalent(nullmod$varCompCov, lmm.genesis$varCompCov)
    expect_equivalent(nullmod$family$family, lmm.genesis$family$family)
    expect_equivalent(nullmod$zeroFLAG, lmm.genesis$zeroFLAG)
    expect_true(all(abs(nullmod$cholSigmaInv - lmm.genesis$cholSigmaInv) < 1e-7))
    expect_equivalent(nullmod$RSS, lmm.genesis$RSS)
})


if (test_binary) {
test_that("fitNullMod matches fitNullMM - binary", {
    ## this never goes FALSE anymore
    ## varCompZero <- TRUE
    ## while(varCompZero){
        dat <- .testNullInputs(binary=TRUE)
        df <- .asDataFrame(dat)
        
        glmm.genesis <- #tryCatch({
            fitNullMM(df, outcome="y", covars = c("X1", "X2", "X3"), covMatList=dat$cor.mat, family="binomial", verbose=FALSE)
        ## }, 
        ## warning = function(w){return(list(message = "warning"))},
        ## error = function(e){return(list(message = "error"))}
        ## )
        ## if (!is.null(glmm.genesis$message)) next
        ## if (glmm.genesis$varComp[1] != 0 ) varCompZero <- FALSE
    ## }

    nullmod <- .fitNullModel(dat$y, dat$X, covMatList=dat$cor.mat, family="binomial", verbose=FALSE)
    
    expect_equivalent(nullmod$fixef, glmm.genesis$fixef, tolerance=1e-6)
    expect_equivalent(nullmod$betaCov, glmm.genesis$betaCov, tolerance=1e-6)
    expect_equivalent(nullmod$resid.marginal, glmm.genesis$resid.marginal, tolerance=1e-6)
    expect_equivalent(nullmod$logLik, glmm.genesis$logLik, tolerance=1e-6)
    expect_equivalent(nullmod$logLikR, glmm.genesis$logLikR, tolerance=1e-6)

    expect_equivalent(nullmod$AIC, glmm.genesis$AIC, tolerance=1e-6)
    expect_equivalent(nullmod$workingY, glmm.genesis$workingY, tolerance=1e-6)
    expect_equivalent(nullmod$model.matrix, glmm.genesis$model.matrix, tolerance=1e-6)
    expect_equivalent(nullmod$varComp, glmm.genesis$varComp, tolerance=1e-6)
    expect_equivalent(nullmod$varCompCov, glmm.genesis$varCompCov, tolerance=1e-4)  ## not smaller then 1e-6 because values were not updated after convergence. However, this is not important, because the variance component values are converged at the desired level!
    expect_equivalent(nullmod$family$family, glmm.genesis$family$family)
    expect_equivalent(nullmod$zeroFLAG, glmm.genesis$zeroFLAG)
    expect_true(all(abs(nullmod$cholSigmaInv - glmm.genesis$cholSigmaInv) < 1e-6))
})
}


test_that("fitNullMod matches fitNullMM - no variance components", {
    n <- 100
    dat <- .testNullInputs(n)
    df <- .asDataFrame(dat)
    
    ## I need to create a junk correlation matrix that will zero out as having variance component zero...  
    varCompJunk <- TRUE
    while(varCompJunk){
	cor.mat <- diag(rnorm(n, sd = 0.01))
	dimnames(cor.mat) <- list(df$scanID, df$scanID)
	lmm.genesis <- fitNullMM(df, "y", covars=c("X1", "X2", "X3"), covMatList=cor.mat, group.var="group", verbose=FALSE)
	if (lmm.genesis$varComp[1] == 0 ) varCompJunk <- FALSE
    }

    nullmod <- .fitNullModel(dat$y, dat$X, group.idx=dat$group.idx, verbose=FALSE)

    expect_equivalent(nullmod$fixef, lmm.genesis$fixef, tolerance=1e-6)
    expect_equivalent(nullmod$betaCov, lmm.genesis$betaCov, tolerance=1e-6)
    expect_equivalent(nullmod$resid.marginal, lmm.genesis$resid.marginal, tolerance=1e-6)
    expect_equivalent(nullmod$logLik, lmm.genesis$logLik, tolerance=1e-6)
    expect_equivalent(nullmod$logLikR, lmm.genesis$logLikR, tolerance=1e-6)

    ## currently GENESIS has a mistake, in AUC calculation it uses the number of 
    ## matrices and groups used, but not the actual number of non-zero variance components. 
    ## so this "2" fixes for it. 
    expect_equivalent(nullmod$AIC, lmm.genesis$AIC - 2, tolerance=1e-6)
    expect_equivalent(nullmod$workingY, lmm.genesis$workingY, tolerance=1e-6)
    expect_equivalent(nullmod$model.matrix, lmm.genesis$model.matrix, tolerance=1e-6)
    expect_equivalent(nullmod$varComp, lmm.genesis$varComp[2:3], tolerance=1e-6)
    expect_equivalent(nullmod$varCompCov, lmm.genesis$varCompCov[2:3, 2:3], tolerance=1e-4)
    expect_equivalent(nullmod$family$family, lmm.genesis$family$family)
    expect_true(all(abs(nullmod$cholSigmaInv - lmm.genesis$cholSigmaInv) < 1e-6))
    expect_equivalent(nullmod$RSS, lmm.genesis$RSS, tolerance=1e-6)
})




test_that("fitNullModel matches fitNulMM", {
    svd <- .testData()
    grm <- .testGRM(svd)
    gen1 <- fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, verbose=FALSE)
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
    expect_true(all(abs(gen2$cholSigmaInv - gen1$cholSigmaInv) < 1e-7))
    expect_equivalent(gen2$RSS, gen1$RSS)

    seqClose(svd)
})

test_that("fitNullModel matches fitNulMM - group", {
    svd <- .testData()
    grm <- .testGRM(svd)
    gen1 <- fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, group.var="status", verbose=FALSE)
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
    expect_true(all(abs(gen2$cholSigmaInv - gen1$cholSigmaInv) < 1e-7))
    expect_equivalent(gen2$RSS, gen1$RSS)

    seqClose(svd)
})

if (test_binary) {
test_that("fitNullModel matches fitNulMM - binary", {
    svd <- .testData()
    grm <- .testGRM(svd)
    gen1 <- fitNullMM(sampleData(svd), outcome="status", covars=c("sex", "age"), covMatList=grm, family="binomial", verbose=FALSE)
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
    lmm.genesis <- fitNullReg(sampleData(svd), outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)

    expect_equivalent(nullmod$fixef, lmm.genesis$fixef)
    expect_equivalent(nullmod$betaCov, lmm.genesis$betaCov)
    expect_equivalent(nullmod$resid.marginal, lmm.genesis$resid.response)
    expect_equivalent(nullmod$logLik, lmm.genesis$logLik)
    expect_equivalent(nullmod$AIC, lmm.genesis$AIC)
    expect_equivalent(nullmod$workingY, lmm.genesis$workingY)
    expect_equivalent(nullmod$model.matrix, lmm.genesis$model.matrix)
    expect_equal(nullmod$family$family, lmm.genesis$family$family)

    seqClose(svd)
})

if (test_binary) {
test_that("fitNullModel matches fitNullReg - binary", {
    svd <- .testData()
    lmm.genesis <- fitNullReg(sampleData(svd), outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
    nullmod <- fitNullModel(svd, outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)

    expect_equivalent(nullmod$fixef, lmm.genesis$fixef)
    expect_equivalent(nullmod$betaCov, lmm.genesis$betaCov)
    expect_equivalent(nullmod$resid.marginal, lmm.genesis$resid.response)
    expect_equivalent(nullmod$logLik, lmm.genesis$logLik)
    expect_equivalent(nullmod$AIC, lmm.genesis$AIC)
    expect_equivalent(nullmod$workingY, lmm.genesis$workingY)
    expect_equivalent(nullmod$model.matrix, lmm.genesis$model.matrix)
    expect_equal(nullmod$family$family, lmm.genesis$family$family)

    seqClose(svd)
})
}


test_that("assocTestSingle matches assocTestMM - Wald", {
    svd <- .testData()
    grm <- .testGRM(svd)
    seqResetFilter(svd, verbose=FALSE)
    nullmod <- fitNullMM(sampleData(svd), outcome="outcome", covMatList=grm, verbose=FALSE)
    
    # multiallelic variants are handled differently
    snv <- isSNV(svd, biallelic=TRUE)
    assoc1 <- assocTestMM(svd, nullmod, test="Wald", snp.include=which(snv), verbose=FALSE)
    
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
    nullmod <- fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, verbose=FALSE)
    
    # multiallelic variants are handled differently
    snv <- isSNV(svd, biallelic=TRUE)
    assoc1 <- assocTestMM(svd, nullmod, test="Score", snp.include=which(snv), verbose=FALSE)
    
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
    nullmod <- fitNullMM(sampleData(svd), outcome="status", covars=c("sex", "age"), covMatList=grm, family="binomial", verbose=FALSE)
    
    # multiallelic variants are handled differently
    snv <- isSNV(svd, biallelic=TRUE)
    assoc1 <- assocTestMM(svd, nullmod, test="Score", snp.include=which(snv), verbose=FALSE)
    
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
    nullmod <- fitNullMM(sampleData(svd), outcome="outcome", covars="sex", covMatList=grm, verbose=FALSE)
    
    # multiallelic variants are handled differently
    assoc1 <- assocTestMM(svd, nullmod, ivars="sex", verbose=FALSE)
    
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



test_that("assocTestSingle matches assocTestMM - GenotypeData", {
    genoData <- .testGenoData()
    grm <- .testGenoDataGRM(genoData)
    chr <- getChromosome(genoData)

    nullmod1 <- fitNullMM(getScanAnnotation(genoData), outcome="outcome", covars="sex", covMatList=grm, verbose=FALSE)
    # can't test Y chr with females
    # also some XY SNPs have exactly the same separation as sex, so can't use it as a covariate
    assoc1 <- assocTestMM(genoData, nullmod1, test="Score", chromosome=c(21:23,26), verbose=FALSE)
    
    iterator <- GenotypeBlockIterator(genoData, snpBlock=1000)
    nullmod2 <- fitNullModel(genoData, outcome="outcome", covars="sex", cov.mat=grm, verbose=FALSE)
    assoc2 <- assocTestSingle(iterator, nullmod2, test="Score", verbose=FALSE)
    assoc2 <- assoc2[assoc2$chr %in% c(21:22, "X", "M"),]

    expect_equal(assoc1$snpID, assoc2$variant.id)
    expect_equal(ifelse(assoc1$minor.allele == "A", assoc1$MAF, 1-assoc1$MAF), assoc2$freq)
    expect_equal(assoc1$Score, assoc2$Score)
    expect_equal(assoc1$Var, (assoc2$Score.SE)^2)
    expect_equal(assoc1$Score.Stat, (assoc2$Score.Stat)^2)
    expect_equal(assoc1$Score.pval, assoc2$Score.pval)

    close(genoData)
})


test_that("assocTestSingle matches assocTestMM - GenotypeData - XY", {
    genoData <- .testGenoData()
    grm <- .testGenoDataGRM(genoData)

    nullmod1 <- fitNullMM(getScanAnnotation(genoData), outcome="outcome", covMatList=grm, verbose=FALSE)
    assoc1 <- assocTestMM(genoData, nullmod1, test="Score", chromosome=24, verbose=FALSE)

    iterator <- GenotypeBlockIterator(genoData, snpBlock=1000)
    nullmod2 <- fitNullModel(genoData, outcome="outcome", cov.mat=grm, verbose=FALSE)
    assoc2 <- assocTestSingle(iterator, nullmod2, verbose=FALSE)
    assoc2 <- assoc2[assoc2$chr == "XY",]

    expect_equal(assoc1$snpID, assoc2$variant.id)
    expect_equal(ifelse(assoc1$minor.allele == "A", assoc1$MAF, 1-assoc1$MAF), assoc2$freq)
    expect_equal(assoc1$Score, assoc2$Score)
    expect_equal(assoc1$Var, (assoc2$Score.SE)^2)
    expect_equal(assoc1$Score.Stat, (assoc2$Score.Stat)^2)
    expect_equal(assoc1$Score.pval, assoc2$Score.pval)

    close(genoData)
})
 

test_that("assocTestSingle matches assocTestMM - GenotypeData - Ychr", {
    genoData <- .testGenoData()
    grm <- .testGenoDataGRM(genoData)

    # select males only
    males <- getScanID(genoData, index=(getSex(genoData) == "M"))
    
    nullmod1 <- fitNullMM(getScanAnnotation(genoData), outcome="outcome", covMatList=grm, scan.include=males, verbose=FALSE)
    assoc1 <- assocTestMM(genoData, nullmod1, test="Score", chromosome=25, verbose=FALSE)

    iterator <- GenotypeBlockIterator(genoData, snpBlock=1000)
    nullmod2 <- fitNullModel(genoData, outcome="outcome", cov.mat=grm, sample.id=males, verbose=FALSE)
    assoc2 <- assocTestSingle(iterator, nullmod2, verbose=FALSE)
    assoc2 <- assoc2[assoc2$chr == "Y",]

    expect_equal(assoc1$snpID, assoc2$variant.id)
    expect_equal(ifelse(assoc1$minor.allele == "A", assoc1$MAF, 1-assoc1$MAF), assoc2$freq)
    expect_equal(assoc1$Score, assoc2$Score)
    expect_equal(assoc1$Var, (assoc2$Score.SE)^2)
    expect_equal(assoc1$Score.Stat, (assoc2$Score.Stat)^2)
    expect_equal(assoc1$Score.pval, assoc2$Score.pval)

    close(genoData)
})
 

test_that(".alleleFreq matches alleleFreq", {
    genoData <- .testGenoData()
    geno <- GWASTools::getGenotype(genoData, transpose=TRUE)
    freq1 <- alleleFreq(geno, chromChar=getChromosome(genoData, char=TRUE), sex=getSex(genoData))
    freq2 <- .alleleFreq(genoData, geno)
    expect_equal(freq1, freq2)
})
   

test_that("admixMap matches admixMapMM", {

# create file with multiple ancestries
gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
tmpfile <- tempfile()
file.copy(gdsfile, tmpfile)
gds <- openfn.gds(tmpfile, readonly=FALSE)
nsnp <- objdesp.gdsn(index.gdsn(gds, "snp.id"))$dim
nsamp <- objdesp.gdsn(index.gdsn(gds, "sample.id"))$dim
dosage_eur <- sample(0:2, nsnp*nsamp, replace=TRUE)
dosage_afr <- ifelse(dosage_eur == 2, 0, sample(0:1, nsnp*nsamp, replace=TRUE))
dosage_amer <- 2 - dosage_eur - dosage_afr
add.gdsn(gds, "dosage_eur", matrix(dosage_eur, nrow=nsamp, ncol=nsnp))
add.gdsn(gds, "dosage_afr", matrix(dosage_afr, nrow=nsamp, ncol=nsnp))
add.gdsn(gds, "dosage_amer", matrix(dosage_amer, nrow=nsamp, ncol=nsnp))
closefn.gds(gds)
        
# read GRM
pcrfile <- system.file("extdata", "HapMap_ASW_MXL_pcrelate.gds", package="GENESIS")
pcr <- openfn.gds(pcrfile)
mypcrel <- pcrelateMakeGRM(pcr)
closefn.gds(pcr)

# generate a phenotype
set.seed(4)
pheno <- rnorm(nsamp, mean = 0, sd = 1)
covar <- sample(0:1, nsamp, replace=TRUE)

# make ScanAnnotationDataFrame
scanAnnot <- GWASTools::ScanAnnotationDataFrame(data.frame(scanID = rownames(mypcrel), 
              covar, pheno, stringsAsFactors=FALSE))

# read in GDS data
gds <- openfn.gds(tmpfile)
genoDataList <- list()
for (anc in c("eur", "afr", "amer")){
  gdsr <- GWASTools::GdsGenotypeReader(gds, genotypeVar=paste0("dosage_", anc))
  genoDataList[[anc]] <- GWASTools::GenotypeData(gdsr, scanAnnot=scanAnnot)
}
     
# fit the null mixed model
nullmod <- fitNullMM(scanData = scanAnnot, outcome = "pheno", covars = "covar", covMatList = mypcrel, verbose=FALSE)

# run the association test
a1 <- admixMapMM(genoDataList, nullMMobj = nullmod, verbose=FALSE)


# new method
nullmod2 <- fitNullModel(scanAnnot, outcome = "pheno", covars = "covar", cov.mat = mypcrel, verbose=FALSE)

# iterators
genoIterators <- lapply(genoDataList, GWASTools::GenotypeBlockIterator)
a2 <- admixMap(genoIterators, nullmod2)

for (v in c("snpID", "chr", "n", "eur.freq", "afr.freq", "amer.freq")) {
    expect_equal(a1[[v]], a2[[v]])
}

## # instability with new null model - same as GxE
## # this worked before we switched to using new calcXtilde and Matrix objects
## nm.tmp <- nullmod
## nm.tmp$sample.id <- nm.tmp$scanID
## lapply(genoIterators, GWASTools::resetIterator)
## a3 <- admixMap(genoIterators, nm.tmp)

## for (v in names(a3)) {
##     expect_equal(a1[[v]], a3[[v]])
## }

GWASTools::close(genoDataList[[1]])
unlink(tmpfile)
})





## defunct functions
    
## test_that("nullModelTestPrep matches calculateProjection", {
##     n <- 100
##     dat <- .testNullInputs(n)
##     geno <- .testGenoMatrix(n)
##     df <- .asDataFrame(dat)
    
##     # basic
##     nullmod <- .fitNullModel(dat$y, dat$X, verbose=FALSE)
##     Xtilde <- calcXtilde(nullmod, geno)
    
##     nullmod.orig <- fitNullReg(df, outcome="y", covars=c("X1", "X2", "X3"), verbose=FALSE)
##     proj <- .calculateProjection(nullmod.orig, test="", burden.test="")

##     expect_true(all(abs(nullmod$Xtilde - crossprod(proj$Mt, geno)) < 1e-7))
##     expect_true(all(abs(nullmod$resid - proj$resid) < 1e-7))
    
##     # with covMatList
##     nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
##     Xtilde <- calcXtilde(nullmod, geno)

##     nullmod.orig <- fitNullMM(df, outcome="y", covars=c("X1", "X2", "X3"), covMatList=dat$cor.mat, verbose=FALSE)
##     proj <- .calculateProjection(nullmod.orig, test="", burden.test="")

##     expect_true(all(abs(Xtilde - crossprod(proj$Mt, geno)) < 1e-7))
##     expect_true(all(abs(nullmod$resid - proj$resid) < 1e-7))
    
##     # with group
##     nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, group.idx=dat$group.idx, verbose=FALSE)
##     Xtilde <- calcXtilde(nullmod, geno)

##     nullmod.orig <- fitNullMM(df, outcome="y", covars=c("X1", "X2", "X3"), covMatList=dat$cor.mat, group.var="group", verbose=FALSE)
##     proj <- .calculateProjection(nullmod.orig, test="", burden.test="")

##     expect_true(all(abs(Xtilde - crossprod(proj$Mt, geno)) < 1e-7))
##     expect_true(all(abs(nullmod$resid - proj$resid) < 1e-7))
## })


## if (test_binary) {
## test_that("nullModelTestPrep vs calculateProjection - binary", {
##     n <- 100
##     geno <- .testGenoMatrix(n)
    
##     ## reps <- 0
##     ## varCompZero <- TRUE
##     ## while(varCompZero & reps < 10){
    
##     dat <- .testNullInputs(n, binary=TRUE)
##     df <- .asDataFrame(dat)	
    
##     ## 	glmm.genesis <- tryCatch({
##     ## 		fitNullMM(df, outcome="y", covars=c("X1", "X2", "X3"), covMatList=cor.mat, family="binomial", verbose=FALSE, maxIter=10)
##     ## 				}, 
##     ## 			warning = function(w){return(list(message = "warning"))},
##     ## 			error = function(e){return(list(message = "error"))}
##     ## 			)
##     ## 	if (!is.null(glmm.genesis$message)) next
##     ## 	if (glmm.genesis$varComp[1] != 0 ) varCompZero <- FALSE
##     ##         reps <- reps + 1
##     ## }
##     ## if (varCompZero) stop("could not generate nonzero varComp")

##     # basic
##     nullmod <- .fitNullModel(dat$y, dat$X, family="binomial", verbose=FALSE)
##     Xtilde <- calcXtilde(nullmod, geno)
    
##     nullmod.orig <- fitNullReg(df, outcome="y", covars=c("X1", "X2", "X3"), family="binomial", verbose=FALSE)
##     proj <- .calculateProjection(nullmod.orig, test="", burden.test="")

##     expect_true(all(abs(Xtilde - crossprod(proj$Mt, geno)) < 1e-7))
##     expect_true(all(abs(nullmod$resid - proj$resid) < 1e-7))
    
##     # with covMatList
##     nullmod <- .fitNullModel(dat$y, dat$X, dat$cor.mat, family="binomial", verbose=FALSE)
##     Xtilde <- calcXtilde(nullmod, geno)

##     nullmod.orig <- fitNullMM(df, outcome="y", covars=c("X1", "X2", "X3"), covMatList=dat$cor.mat, family="binomial", verbose=FALSE)
##     proj <- .calculateProjection(nullmod.orig, test="", burden.test="")

##     expect_true(all(abs(Xtilde - crossprod(proj$Mt, geno)) < 1e-7))
##     expect_true(all(abs(nullmod$resid - proj$resid) < 1e-7))
## })
## }
    
## test_that("assocTestAggregate matches assocTestSeqWindow - Burden, Wald", {
##     svd <- .testData()
##     nullmod <- fitNullReg(sampleData(svd), outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
##     assoc1 <- assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="Burden", burden.test="Wald", verbose=FALSE, chromosome=1)
    
##     nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
##     seqSetFilterChrom(svd, include=1, verbose=FALSE)
##     iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
##     assoc2 <- assocTestAggregate(iterator, nullmod, test="Burden", burden.test="Wald", verbose=FALSE)

##     var1 <- assoc1$variantInfo[!(assoc1$variantInfo$freq %in% c(0,1)),]
##     var2 <- do.call(rbind, assoc2$variantInfo)
##     var2 <- var2[!duplicated(paste(var2$variant.id, var2$allele.index)),]
##     expect_equal(nrow(var1), nrow(var2))  
##     expect_equal(var1$variantID, var2$variant.id)
##     expect_equal(var1$allele, var2$allele.index)
##     expect_equal(var1$chr, var2$chr)
##     expect_equal(var1$pos, var2$pos)
##     expect_equal(var1$n.obs, var2$n.obs)
##     expect_equal(var1$freq, var2$freq)
##     expect_equal(var1$weight, var2$weight)

##     res1 <- assoc1$results[assoc1$results$dup %in% 0,]
##     res2 <- assoc2$results[assoc2$results$n.site > 0,]
##     expect_equal(res1$n.site, res2$n.site)
##     expect_equal(res1$Est, res2$Est)
##     expect_equal(res1$SE, res2$Est.SE)
##     expect_equal(res1$Wald.stat, (res2$Wald.Stat)^2)
##     expect_equal(res1$Wald.pval, res2$Wald.pval)
    
##     seqClose(svd)
## })


## test_that("assocTestAggregate matches assocTestSeqWindow - Burden, Score", {
##     svd <- .testData()
##     nullmod <- fitNullReg(sampleData(svd), outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
##     assoc1 <- assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="Burden", burden.test="Score", verbose=FALSE, chromosome=1)
    
##     nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
##     seqSetFilterChrom(svd, include=1, verbose=FALSE)
##     iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
##     assoc2 <- assocTestAggregate(iterator, nullmod, test="Burden", burden.test="Score", verbose=FALSE)

##     res1 <- assoc1$results[assoc1$results$dup %in% 0,]
##     res2 <- assoc2$results[assoc2$results$n.site > 0,]
##     expect_equal(res1$n.site, res2$n.site)
##     expect_equal(res1$Score, res2$Score)
##     expect_equal(res1$Var, (res2$Score.SE)^2)
##     expect_equal(res1$Score.stat, (res2$Score.Stat)^2)
##     expect_equal(res1$Score.pval, res2$Score.pval)
    
##     seqClose(svd)
## })


## test_that("assocTestAggregate matches assocTestSeqWindow - Burden, Wald, LMM", {
##     svd <- .testData()
##     grm <- .testGRM(svd)
##     nullmod <- fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, verbose=FALSE)
##     assoc1 <- assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="Burden", burden.test="Wald", verbose=FALSE, chromosome=1)
    
##     nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), cov.mat=grm, verbose=FALSE)
##     seqSetFilterChrom(svd, include=1, verbose=FALSE)
##     iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
##     assoc2 <- assocTestAggregate(iterator, nullmod, test="Burden", burden.test="Wald", verbose=FALSE)

##     res1 <- assoc1$results[assoc1$results$dup %in% 0,]
##     res2 <- assoc2$results[assoc2$results$n.site > 0,]
##     expect_equal(res1$n.site, res2$n.site)
##     expect_equal(res1$Est, res2$Est)
##     expect_equal(res1$SE, res2$Est.SE)
##     expect_equal(res1$Wald.stat, (res2$Wald.Stat)^2)
##     expect_equal(res1$Wald.pval, res2$Wald.pval)
    
##     seqClose(svd)
## })


## test_that("assocTestAggregate matches assocTestSeqWindow - Burden, Score, LMM", {
##     svd <- .testData()
##     grm <- .testGRM(svd)
##     nullmod <- fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, verbose=FALSE)
##     assoc1 <- assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="Burden", burden.test="Score", verbose=FALSE, chromosome=1)
    
##     nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), cov.mat=grm, verbose=FALSE)
##     seqSetFilterChrom(svd, include=1, verbose=FALSE)
##     iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
##     assoc2 <- assocTestAggregate(iterator, nullmod, test="Burden", burden.test="Score", verbose=FALSE)

##     res1 <- assoc1$results[assoc1$results$dup %in% 0,]
##     res2 <- assoc2$results[assoc2$results$n.site > 0,]
##     expect_equal(res1$n.site, res2$n.site)
##     expect_equal(res1$Score, res2$Score)
##     expect_equal(res1$Var, (res2$Score.SE)^2)
##     expect_equal(res1$Score.stat, (res2$Score.Stat)^2)
##     expect_equal(res1$Score.pval, res2$Score.pval)
    
##     seqClose(svd)
## })

## test_that("assocTestAggregate matches assocTestSeqWindow - Burden, Binary", {
##     svd <- .testData()
##     nullmod <- fitNullReg(sampleData(svd), outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
##     assoc1 <- assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="Burden", burden.test="Score", verbose=FALSE, chromosome=1)
    
##     nullmod <- fitNullModel(svd, outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
##     seqSetFilterChrom(svd, include=1, verbose=FALSE)
##     iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
##     assoc2 <- assocTestAggregate(iterator, nullmod, test="Burden", burden.test="Score", verbose=FALSE)

##     res1 <- assoc1$results[assoc1$results$dup %in% 0,]
##     res2 <- assoc2$results[assoc2$results$n.site > 0,]
##     expect_equal(res1$n.site, res2$n.site)
##     expect_equal(res1$Score, res2$Score, tolerance=1e-6)
##     expect_equal(res1$Var, (res2$Score.SE)^2, tolerance=1e-6)
##     expect_equal(res1$Score.stat, (res2$Score.Stat)^2, tolerance=1e-6)
##     expect_equal(res1$Score.pval, res2$Score.pval, tolerance=1e-6)
    
##     seqClose(svd)
## })


## if (test_binary) {
## test_that("assocTestAggregate matches assocTestSeqWindow - Burden, Binary, LMM", {
##     svd <- .testData()
##     grm <- .testGRM(svd)
##     nullmod <- fitNullMM(sampleData(svd), outcome="status", covars=c("sex", "age"), covMatList=grm, family="binomial", verbose=FALSE)
##     assoc1 <- assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="Burden", burden.test="Score", verbose=FALSE, chromosome=1)
    
##     nullmod <- fitNullModel(svd, outcome="status", covars=c("sex", "age"), cov.mat=grm, family="binomial", verbose=FALSE)
##     seqSetFilterChrom(svd, include=1, verbose=FALSE)
##     iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
##     assoc2 <- assocTestAggregate(iterator, nullmod, test="Burden", burden.test="Score", verbose=FALSE)

##     res1 <- assoc1$results[assoc1$results$dup %in% 0,]
##     res2 <- assoc2$results[assoc2$results$n.site > 0,]
##     expect_equal(res1$n.site, res2$n.site)
##     expect_equal(res1$Score, res2$Score, tolerance=1e-6)
##     expect_equal(res1$Var, (res2$Score.SE)^2, tolerance=1e-6)
##     expect_equal(res1$Score.stat, (res2$Score.Stat)^2, tolerance=1e-6)
##     expect_equal(res1$Score.pval, res2$Score.pval, tolerance=1e-6)
    
##     seqClose(svd)
## })
## }


## test_that("assocTestAggregate matches assocTestSeqWindow - SKAT", {
##     svd <- .testData()
##     nullmod <- fitNullReg(sampleData(svd), outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
##     assoc1 <- assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="SKAT", verbose=FALSE, chromosome=1)
    
##     nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
##     seqSetFilterChrom(svd, include=1, verbose=FALSE)
##     iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
##     assoc2 <- assocTestAggregate(iterator, nullmod, test="SKAT", verbose=FALSE)

##     res1 <- assoc1$results[assoc1$results$dup %in% 0,]
##     res2 <- assoc2$results[assoc2$results$n.site > 0,]
##     expect_equal(res1$n.site, res2$n.site)
##     expect_equal(res1$Q_0, res2$Q_0)
##     expect_equal(res1$pval_0, res2$pval_0)
##     expect_equal(res1$err_0, res2$err_0)
    
##     seqClose(svd)
## })


## if (test_binary) {
## test_that("assocTestAggregate matches assocTestSeqWindow - SKAT, binary", {
##     svd <- .testData()
##     nullmod <- fitNullReg(sampleData(svd), outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
##     assoc1 <- assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="SKAT", verbose=FALSE, chromosome=1)
    
##     nullmod <- fitNullModel(svd, outcome="status", covars=c("sex", "age"), family="binomial", verbose=FALSE)
##     seqSetFilterChrom(svd, include=1, verbose=FALSE)
##     iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
##     assoc2 <- assocTestAggregate(iterator, nullmod, test="SKAT", verbose=FALSE)

##     res1 <- assoc1$results[assoc1$results$dup %in% 0,]
##     res2 <- assoc2$results[assoc2$results$n.site > 0,]
##     expect_equal(res1$n.site, res2$n.site)
##     expect_equal(res1$Q_0, res2$Q_0, tolerance=1e-6)
##     expect_equal(res1$pval_0, res2$pval_0, tolerance=1e-6)
##     expect_equal(res1$err_0, res2$err_0, tolerance=1e-6)
    
##     seqClose(svd)
## })
## }


## test_that("assocTestAggregate matches assocTestSeqWindow - SKAT, LMM", {
##     svd <- .testData()
##     grm <- .testGRM(svd)
##     nullmod <- fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, verbose=FALSE)
##     assoc1 <- assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="SKAT", verbose=FALSE, chromosome=1)
    
##     nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), cov.mat=grm, verbose=FALSE)
##     seqSetFilterChrom(svd, include=1, verbose=FALSE)
##     iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
##     assoc2 <- assocTestAggregate(iterator, nullmod, test="SKAT", verbose=FALSE)

##     res1 <- assoc1$results[assoc1$results$dup %in% 0,]
##     res2 <- assoc2$results[assoc2$results$n.site > 0,]
##     expect_equal(res1$n.site, res2$n.site)
##     expect_equal(res1$Q_0, res2$Q_0)
##     expect_equal(res1$pval_0, res2$pval_0)
##     expect_equal(res1$err_0, res2$err_0)
    
##     seqClose(svd)
## })


## if (test_binary) {
## test_that("assocTestAggregate matches assocTestSeqWindow - SKAT, binary, LMM", {
##     svd <- .testData()
##     grm <- .testGRM(svd)
##     nullmod <- fitNullMM(sampleData(svd), outcome="status", covars=c("sex", "age"), covMatList=grm, family="binomial", verbose=FALSE)
##     assoc1 <- assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="SKAT", verbose=FALSE, chromosome=1)
    
##     nullmod <- fitNullModel(svd, outcome="status", covars=c("sex", "age"), cov.mat=grm, family="binomial", verbose=FALSE)
##     seqSetFilterChrom(svd, include=1, verbose=FALSE)
##     iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
##     assoc2 <- assocTestAggregate(iterator, nullmod, test="SKAT", verbose=FALSE)

##     res1 <- assoc1$results[assoc1$results$dup %in% 0,]
##     res2 <- assoc2$results[assoc2$results$n.site > 0,]
##     expect_equal(res1$n.site, res2$n.site)
##     expect_equal(res1$Q_0, res2$Q_0, tolerance=1e-6)
##     expect_equal(res1$pval_0, res2$pval_0, tolerance=1e-6)
##     expect_equal(res1$err_0, res2$err_0, tolerance=1e-6)
    
##     seqClose(svd)
## })
## }


## test_that("assocTestAggregate matches assocTestSeqWindow - SKAT-O, LMM", {
##     svd <- .testData()
##     grm <- .testGRM(svd)
##     nullmod <- fitNullMM(sampleData(svd), outcome="outcome", covars=c("sex", "age"), covMatList=grm, verbose=FALSE)
##     assoc1 <- assocTestSeqWindow(svd, nullmod, window.size=100, window.shift=50, test="SKAT", rho=c(0,0.5,1), verbose=FALSE, chromosome=1)
    
##     nullmod <- fitNullModel(svd, outcome="outcome", covars=c("sex", "age"), cov.mat=grm, verbose=FALSE)
##     seqSetFilterChrom(svd, include=1, verbose=FALSE)
##     iterator <- SeqVarWindowIterator(svd, windowSize=1e5, windowShift=5e4, verbose=FALSE)
##     assoc2 <- assocTestAggregate(iterator, nullmod, test="SKAT", rho=c(0,0.5,1), verbose=FALSE)

##     res1 <- assoc1$results[assoc1$results$dup %in% 0,]
##     res2 <- assoc2$results[assoc2$results$n.site > 0,]
##     expect_equal(res1$n.site, res2$n.site)
##     expect_equal(res1$Q_0, res2$Q_0)
##     expect_equal(res1$pval_0, res2$pval_0)
##     expect_equal(res1$err_0, res2$err_0)
    
##     seqClose(svd)
## })

