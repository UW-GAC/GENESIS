context("check single variant association tests")

test_that("singleVarTest - linear, with group", {
    n <- 100
    dat <- .testNullInputs(n)
    geno <- .testGenoMatrix(n)
    
    nullmod <- .fitNullModel(dat$y, dat$X, group.idx=dat$group.idx, verbose=FALSE)
    test1 <- testGenoSingleVar(nullmod, G = geno)
    
    # compare to weighted least squares using weighted lm:	
    res.lm <- numeric(nrow(test1))
    for (i in 1:ncol(geno)){
        lm.temp <- lm(dat$y ~ -1 + dat$X + geno[,i], weights = c(rep(1/nullmod$varComp[1], n/2), 1/rep(nullmod$varComp[2], n/2)))
        res.lm[i] <- summary(lm.temp)$coef[4,"Estimate"]
    }
    
    expect_equal(res.lm, test1$Est)
})


test_that("singleVarTest - linear, without group", {
    n <- 100
    dat <- .testNullInputs(n)
    geno <- .testGenoMatrix(n)
    
    # without group
    nullmod <- .fitNullModel(dat$y, dat$X, verbose=FALSE)
    test1 <- testGenoSingleVar(nullmod, G = geno)
    
    res.lm <- numeric(nrow(test1))
    for (i in 1:ncol(geno)){
        lm.temp <- lm(dat$y ~ -1 + dat$X + geno[,i], )
        res.lm[i] <- summary(lm.temp)$coef[4,"Estimate"]
    }
    
    expect_equal(res.lm, test1$Est)
})


test_that("singleVarTest - logistic - wald", {
    n <- 100
    dat <- .testNullInputs(n, binary=TRUE)
    geno <- .testGenoMatrix(n)
    
    ##comparing the Wald test - in genesis computed using the fixed effects of the null model. 
    test.wald <- data.frame(Est = rep(NA, ncol(geno)), Est.SE = NA, Wald.Stat = NA, Wald.pval = NA)
    res.glm <- test.wald	
    
    for (i in 1:ncol(geno)){
        nullmod <- .fitNullModel(dat$y, cbind(dat$X, geno[,i]), family=binomial(), verbose=FALSE)
        test.wald[i,] <- nullmod$fixef[4,]
        glm.temp <-  glm(dat$y ~ -1 + dat$X + geno[,i], family=binomial())
        res.glm[i,] <- summary(glm.temp)$coef[4,]
    }

    expect_equal(res.glm, test.wald)
    
    ## check that we get appropriate error when using the wald test instead of score with binomial outcomes:
    nullmod <- .fitNullModel(dat$y, dat$X, family=binomial(), verbose=FALSE)
    test.score <- testGenoSingleVar(nullmod, G = geno)

    expect_equal(test.score$Score.pval, test.wald$Wald.pval, tolerance = 0.01)
})


test_that("singleVarTest - logistic - score", {
    n <- 100
    dat <- .testNullInputs(n, binary=TRUE)
    geno <- .testGenoMatrix(n)
    
    nullmod <- .fitNullModel(dat$y, dat$X, family=binomial(), verbose=FALSE)
    test.score <- testGenoSingleVar(nullmod, G = geno, test = "Score")

    res.glm <- rep(NA, ncol(geno))	
    for (i in 1:ncol(geno)){
        g <- geno[,i]
        glm.temp <-  glm(dat$y ~ -1 + dat$X + g, family=binomial())
        rao <- anova(glm.temp, test="Rao")
        res.glm[i] <- rao["g", "Pr(>Chi)"]
    }
    
    expect_equivalent(res.glm, test.score$Score.pval, tolerance=1e-4)
})


### BR and CMP unit tests -- no equivalent built-in functions available
test_that("singleVarTest - logistic - BR", {
    n <- 100
    nullmod <- .testNullmod(n, MM=FALSE, binary=TRUE)
    geno <- .testGenoMatrix(n)
    test.br <- testGenoSingleVar(nullmod, G = geno, test = "BinomiRare")
    expect_true(all(colSums(geno) >= test.br$n.carrier))
    expect_true(all(test.br$n.carrier >= test.br$n.D.carrier))
    expect_true(all(test.br$pval > test.br$mid.pval))

    test.score <- testGenoSingleVar(nullmod, G = geno, test = "Score")
    test.br <- testGenoSingleVar(nullmod, G = geno, test = "BinomiRare", recalc.pval.thresh = 0)
    expect_equal(test.score$Score.pval, test.br$pval)
    expect_equal(test.score$Score.pval, test.br$mid.pval)
})

test_that("singleVarTest - logistic - CMP", {
    n <- 100
    nullmod <- .testNullmod(n, MM=TRUE, binary=TRUE)
    geno <- .testGenoMatrix(n)
    test.cmp <- testGenoSingleVar(nullmod, G = geno, test = "CMP")
    expect_true(all(colSums(geno) >= test.cmp$n.carrier))
    expect_true(all(test.cmp$n.carrier >= test.cmp$n.D.carrier))
    expect_true(all(test.cmp$pval > test.cmp$mid.pval))

    test.score <- testGenoSingleVar(nullmod, G = geno, test = "Score")
    test.cmp <- testGenoSingleVar(nullmod, G = geno, test = "CMP", recalc.pval.thresh = 0)
    expect_equal(test.score$Score.pval, test.cmp$pval)
    expect_equal(test.score$Score.pval, test.cmp$mid.pval)

    # CMP switches to BinomiRare if not a mixed model
    nullmod <- .testNullmod(n, MM=FALSE, binary=TRUE)
    test.cmp <- testGenoSingleVar(nullmod, G = geno, test = "CMP")
    test.br <- testGenoSingleVar(nullmod, G = geno, test = "BinomiRare")
    expect_equal(test.br$pval, test.cmp$pval)
})


test_that("GxE", {
    n <- 100
    dat <- .testNullInputs(n)
    colnames(dat$X) <- c("a", "b", "c")
    geno <- .testGenoMatrix(n)
    
    nullmod <- .fitNullModel(dat$y, dat$X, verbose=FALSE)
    test.gxe <- testGenoSingleVar(nullmod, G = geno, E = dat$X[,3,drop=FALSE], GxE.return.cov = TRUE)
    expect_true(all(c("Est.G:c", "SE.G:c") %in% names(test.gxe$res)))
    expect_equal(length(test.gxe$GxEcovMatList), ncol(geno))

    res.lm <- test.gxe$res[,1:4]
    tmp <- .asDataFrame(dat)
    for (i in 1:ncol(geno)){
        lm.temp <- lm("y ~ b + c + g + c:g", data=cbind(tmp, g=geno[,i]))
        res.lm[i,"Est.G"] <- summary(lm.temp)$coef["g",1]
        res.lm[i,"SE.G"] <- summary(lm.temp)$coef["g",2]
        res.lm[i,"Est.G:c"] <- summary(lm.temp)$coef["c:g",1]
        res.lm[i,"SE.G:c"] <- summary(lm.temp)$coef["c:g",2]
    }
    
    expect_equal(res.lm$Est.G, test.gxe$res$Est.G)
    expect_equal(res.lm$SE.G, test.gxe$res$SE.G)
    expect_equal(res.lm$`Est.G:c`, test.gxe$res$`Est.G:c`)
    expect_equal(res.lm$`SE.G:c`, test.gxe$res$`SE.G:c`)
})

test_that("singleVarTest - SPA", {
    n <- 100
    dat <- .testNullInputs(n, binary=TRUE)
    geno <- .testGenoMatrix(n)
    
    nullmod <- .fitNullModel(dat$y, dat$X, family=binomial(), verbose=FALSE)
    test.score <- testGenoSingleVar(nullmod, G = geno, test = "Score")
    test.spa <- testGenoSingleVar(nullmod, G = geno, test = "Score.SPA")
    expect_true(max(test.score$Score.pval - test.spa$SPA.pval) < 0.01)
})

test_that("SPA_pval works with empty input", {
    n <- 100
    dat <- .testNullInputs(n, binary=TRUE)
    geno <- .testGenoMatrix(n)
    
    nullmod <- .fitNullModel(dat$y, dat$X, family=binomial(), verbose=FALSE)
    empty <- data.frame(Score=numeric(), Score.SE=numeric(), Score.Stat=numeric(), Score.pval=numeric())
    empty.G <- matrix(nrow=n, ncol=0)
    test.spa <- SPA_pval(score.result=empty, nullmod=nullmod, G=empty.G)
    expect_equal(nrow(test.spa), 0)
})


test_that("small null model", {
    n <- 100
    dat <- .testNullInputs(n, binary=TRUE)
    nullmod.big <- .fitNullModel(dat$y, dat$X, covMatList=dat$cor.mat, family=binomial(),
                                 return.small=FALSE, verbose=FALSE)
    geno <- .testGenoMatrix(n)
    test.br.big <- testGenoSingleVar(nullmod.big, G = geno, test = "BinomiRare")
    test.cmp.big <- testGenoSingleVar(nullmod.big, G = geno, test = "CMP")
    
    nullmod.small1 <- .fitNullModel(dat$y, dat$X, covMatList=dat$cor.mat, family=binomial(),
                                    return.small=TRUE, verbose=FALSE)
    expect_equal(setdiff(names(nullmod.big), names(nullmod.small1)),
                 c("cholSigmaInv", "CX", "CXCXI"))
    
    test.br.small1 <- testGenoSingleVar(nullmod.small1, G = geno, test = "BinomiRare")
    expect_equal(test.br.big, test.br.small1)
    test.cmp.small1 <- testGenoSingleVar(nullmod.small1, G = geno, test = "CMP")
    expect_equal(test.cmp.big, test.cmp.small1)

    nullmod.small2 <- nullModelSmall(nullmod.big)
    expect_equal(setdiff(names(nullmod.big), names(nullmod.small2)),
                 c("cholSigmaInv", "CX", "CXCXI"))

    test.br.small2 <- testGenoSingleVar(nullmod.small2, G = geno, test = "BinomiRare")
    expect_equal(test.br.big, test.br.small2)
    test.cmp.small2 <- testGenoSingleVar(nullmod.small2, G = geno, test = "CMP")
    expect_equal(test.cmp.big, test.cmp.small2)

    # expect_error(testGenoSingleVar(nullmod.small2, G = geno, test = "BinomiRare",
    #                                recalc.pval.thresh=0.5),
    #              "small null.model cannot be used with test options provided")
    # expect_error(testGenoSingleVar(nullmod.small2, G = geno, test = "Score"),
    #              "small null.model cannot be used with test options provided")

    expect_error(.checkNullModelTestSingle(nullmod.small2, test = "BinomiRare", 
        recalc.pval.thresh = 0.5, fast.score.SE = FALSE, GxE = NULL), 
        "small null.model cannot be used with test options provided")
    expect_error(.checkNullModelTestSingle(nullmod.small2, test = "Score", 
        recalc.pval.thresh = 1, fast.score.SE = FALSE, GxE = NULL), 
        "small null.model cannot be used with test options provided")
    expect_error(.checkNullModelTestSingle(nullmod.small2, test = "Score.SPA", 
        recalc.pval.thresh = 1, fast.score.SE = FALSE, GxE = NULL), 
        "small null.model cannot be used with test options provided")
    expect_error(.checkNullModelTestSingle(nullmod.small2, test = "Score", 
        recalc.pval.thresh = 1, fast.score.SE = TRUE, GxE = NULL), 
        "null.model must have se.correction when fast.score.SE = TRUE; re-fit your null.model using `fitNullModelFastScore` or update your null.model using `nullModelFastScore`")
})
