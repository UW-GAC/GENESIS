context("check single variant association tests")

test_that("singleVarTest - linear, with group", {
    n <- 100
    dat <- .testNullInputs(n)
    geno <- .testGenoMatrix(n)
    
    nullmod <- .fitNullModel(dat$y, dat$X, group.idx=dat$group.idx, verbose=FALSE)
    test.wald <- testGenoSingleVar(nullmod, G = geno, test = "Wald")
    
    # compare to weighted least squares using weighted lm:	
    res.lm <- test.wald
    for (i in 1:ncol(geno)){
        lm.temp <- lm(dat$y ~ -1 + dat$X + geno[,i], weights = c(rep(1/nullmod$varComp[1], n/2), 1/rep(nullmod$varComp[2], n/2)))
        res.lm[i,] <- summary(lm.temp)$coef[4,]
    }
    
    expect_equal(res.lm$Est, test.wald$Est)
    expect_equal(res.lm$Wald.Stat, test.wald$Wald.Stat)
})


test_that("singleVarTest - linear, without group", {
    n <- 100
    dat <- .testNullInputs(n)
    geno <- .testGenoMatrix(n)
    
    # without group
    nullmod <- .fitNullModel(dat$y, dat$X, verbose=FALSE)
    test.wald <- testGenoSingleVar(nullmod, G = geno, test = "Wald")
    
    res.lm <- test.wald
    for (i in 1:ncol(geno)){
        lm.temp <- lm(dat$y ~ -1 + dat$X + geno[,i], )
        res.lm[i,] <- summary(lm.temp)$coef[4,]
    }
    
    expect_equal(res.lm$Est, test.wald$Est)
    expect_equal(res.lm$Wald.Stat, test.wald$Wald.Stat)
})


test_that("singleVarTest - logistic - wald", {
    n <- 100
    dat <- .testNullInputs(n, binary=TRUE)
    geno <- .testGenoMatrix(n)
    
    ##comparing the Wald test - in genesis computed using the fixed effects of the null model. 
    test.wald <- data.frame(Est = rep(NA, ncol(geno)), Est.SE = NA, Wald.Stat = NA, Wald.pval = NA)
    res.glm <- test.wald	
    
    for (i in 1:ncol(geno)){
        nullmod <- .fitNullModel(dat$y, cbind(dat$X, geno[,i]), family="binomial", verbose=FALSE)
        test.wald[i,] <- nullmod$fixef[4,]
        glm.temp <-  glm(dat$y ~ -1 + dat$X + geno[,i], family="binomial")
        res.glm[i,] <- summary(glm.temp)$coef[4,]
    }

    expect_equal(res.glm, test.wald)
    
    ## check that we get appropriate error when using the wald test instead of score with binomial outcomes:
    nullmod <- .fitNullModel(dat$y, dat$X, family="binomial", verbose=FALSE)
    expect_message(test.score <- testGenoSingleVar(nullmod, G = geno, test = "Wald"), "Cannot use Wald test")

    expect_equal(colnames(test.score)[1], "Score")
    expect_equal(test.score$Score.pval, test.wald$Wald.pval, tolerance = 0.01)
})


test_that("singleVarTest - logistic - score", {
    n <- 100
    dat <- .testNullInputs(n, binary=TRUE)
    geno <- .testGenoMatrix(n)
    
    nullmod <- .fitNullModel(dat$y, dat$X, family="binomial", verbose=FALSE)
    test.score <- testGenoSingleVar(nullmod, G = geno, test = "Score")

    res.glm <- rep(NA, ncol(geno))	
    for (i in 1:ncol(geno)){
        g <- geno[,i]
        glm.temp <-  glm(dat$y ~ -1 + dat$X + g, family="binomial")
        rao <- anova(glm.temp, test="Rao")
        res.glm[i] <- rao["g", "Pr(>Chi)"]
    }
    
    expect_equivalent(res.glm, test.score$Score.pval, tolerance=1e-4)
})


test_that("GxE", {
    n <- 100
    dat <- .testNullInputs(n)
    colnames(dat$X) <- c("a", "b", "c")
    geno <- .testGenoMatrix(n)
    
    nullmod <- .fitNullModel(dat$y, dat$X, verbose=FALSE)
    test.gxe <- testGenoSingleVar(nullmod, G = geno, E = dat$X[,3,drop=FALSE], test = "Wald", GxE.return.cov = TRUE)
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

    expect_message(test.gxe <- testGenoSingleVar(nullmod, G = geno, E = dat$X[,3,drop=FALSE], test = "Score"))
})
