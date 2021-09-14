.iterateDGLMM <- function(){

  if (is.null(group.idx) group.idx <- list(resid.var = 1:length(y))
  # fit the initial LMM assuming no dcovars effects
  # all of the dispersion parameters (dmu; i.e. phi) = 1
  dmu <- rep(1, length(y))
  newstart <- start
  LMMreps <- 0

  repeat({
    LMMreps <- LMMreps + 1

    # fit the LMM
    vc.mod <- .runAIREMLgaussian(y, X, start = newstart, covMatList = covMatList,
                                 group.idx = group.idx, diagV = dmu,
                                 AIREML.tol = AIREML.tol, drop.zeros = drop.zeros,
                                 max.iter = max.iter, EM.iter = EM.iter, verbose = verbose)

    # calculate the deviance residuals
    mu <- family$linkinv(vc.mod$eta) # = eta for gaussian
    dev <- family$dev.resids(y, mu) # = (y - mu)^2 for gaussian

    # regress the deviance values on the dcovars
    ### option 1 - use glm
    d.mod <- glm(dev ~ dX, family = dfamily)

    deta <- d.mod$linear.predictors
    dmu <- d.mod$fitted.values
    ###

    ### option 2 - create the working vector and the
    deta <- dfamily$linkfun(dmu)
    dvmu <- dfamily$variance(dmu) # = dmu^2 for Gamma
    dgmuinv <- dfamily$mu.eta(deta)

    workingDev <- deta + (dev - dmu)/dgmuinv
    dW <- dgmuinv^2/dvmu
    d.fit <- lm.wfit(dX, workingDev, dW)

    deta <- d.fit$fitted.values
    dmu <- dfamily$linkinv(deta)
    ###

    # check for convergence
    if(LMMreps > 1){
      if(sqrt(sum((vc.mod$eta - eta)^2)) < AIREML.tol && sqrt(sum((d.mod$linear.predictors - deta)^2)) < AIREML.tol){
         converged <- TRUE
         (break)()
      }
    }

    # check if exceeded the number of iterations
    if(LMMreps == max.iter){
      vc.mod$converged <- FALSE
      warning("Maximum number of iterations for DGLM reached without convergence!")
      (break)()
    }

    # update starting values to current variance component estimates
    newstart <- vc.mod$varComp
    newstart[vc.mod$zeroFLAG] <- AIREML.tol
    # update eta
    eta <- vc.mod$eta
    # update deta
    deta <-

  })




  vc.mod <- .runAIREMLgaussian(y, X, start = start, covMatList = covMatList,
                               group.idx = group.idx, diagV = diagV,
                               AIREML.tol = AIREML.tol, drop.zeros = drop.zeros,
                               max.iter = max.iter, EM.iter = EM.iter, verbose = verbose)

}
