context("check update of null model after ranknormalizing and re-scaling residuals")

test_that("updateOutcome", {
    dat <- .testNullInputs()
    group.idx <- dat$group.idx
    cor.mat <- dat$cor.mat
    covMatList <- list(A = cor.mat)

    nullmod <- .fitNullModel(dat$y, dat$X, group.idx = group.idx, covMatList, verbose=FALSE)

    group.ind <- 1
    expect_equal(.averageGroupVar(nullmod$varComp, covMatList = covMatList, group.idx = group.idx[group.ind]),
                 nullmod$varComp[1]*mean(diag(cor.mat)[group.idx[[group.ind]]]) + nullmod$varComp[group.ind + 1])

    group.ind <- 2
    expect_equal(.averageGroupVar(nullmod$varComp, covMatList = covMatList, group.idx = group.idx[group.ind]),
                 nullmod$varComp[1]*mean(diag(cor.mat)[group.idx[[group.ind]]]) + nullmod$varComp[group.ind + 1])

    expect_equal(.averageGroupVar(nullmod$varComp, covMatList = NULL, group.idx = group.idx[group.ind]),
                 nullmod$varComp[group.ind + 1])
    
    throws_error(.averageGroupVar(nullmod$varComp, covMatList = NULL, group.idx = group.idx))
    throws_error(.averageGroupVar(nullmod$varComp, covMatList = NULL, group.idx = NULL))
    throws_error(.averageGroupVar(nullmod$varComp, covMatList = NULL, group.idx = group.idx[[1]]))

    nullmod2 <- updateNullModOutcome(nullmod, covMatList = covMatList, rankNorm.option = c("by.group"), rescale = c("none"), verbose=FALSE)


    expect_true(abs(.averageGroupVar(nullmod2$varComp, covMatList = covMatList, group.idx = group.idx[1]) - 1) < 0.1 )
    expect_true(abs(.averageGroupVar(nullmod2$varComp, covMatList = covMatList, group.idx = group.idx[2]) - 1) < 0.1 )

    # some comparisons when the residuals are rank-normalized and re-scaled within gorups
    average.group1.model.var.nullmod <- .averageGroupVar(nullmod$varComp, covMatList = covMatList, group.idx = group.idx[1])
    average.group2.model.var.nullmod <- .averageGroupVar(nullmod$varComp, covMatList = covMatList, group.idx = group.idx[2])

    nullmod3 <- updateNullModOutcome(nullmod, covMatList = covMatList, rankNorm.option = c("by.group"), rescale = c("residSD"), verbose=FALSE)
    nullmod4 <- updateNullModOutcome(nullmod, covMatList = covMatList, rankNorm.option = c("by.group"), rescale = c("model"), verbose=FALSE)

    average.group1.model.var.nullmod3 <- .averageGroupVar(nullmod3$varComp, covMatList = covMatList, group.idx = group.idx[1])
    average.group2.model.var.nullmod3 <- .averageGroupVar(nullmod3$varComp, covMatList = covMatList, group.idx = group.idx[2])
    average.group1.model.var.nullmod4 <- .averageGroupVar(nullmod4$varComp, covMatList = covMatList, group.idx = group.idx[1])
    average.group2.model.var.nullmod4 <- .averageGroupVar(nullmod4$varComp, covMatList = covMatList, group.idx = group.idx[2])

    # compare to model 3:
    expect_true(average.group1.model.var.nullmod3/average.group1.model.var.nullmod < 1.1 & average.group1.model.var.nullmod3/average.group1.model.var.nullmod > 0.9 )
    expect_true(average.group2.model.var.nullmod3/average.group2.model.var.nullmod < 1.1 & average.group2.model.var.nullmod3/average.group2.model.var.nullmod > 0.9 )

    # compare to model 4:
    expect_true(average.group1.model.var.nullmod4/average.group1.model.var.nullmod < 1.1 & average.group1.model.var.nullmod4/average.group1.model.var.nullmod > 0.9 )
    expect_true(average.group2.model.var.nullmod4/average.group2.model.var.nullmod < 1.1 & average.group2.model.var.nullmod4/average.group2.model.var.nullmod > 0.9 )
})
