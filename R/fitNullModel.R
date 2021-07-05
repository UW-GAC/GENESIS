
setGeneric("fitNullModel", function(x, ...) standardGeneric("fitNullModel"))

setMethod("fitNullModel",
          "data.frame",
          function(x, outcome,
                   covars = NULL,
                   cov.mat = NULL,
                   group.var = NULL,
                   family = "gaussian",
                   two.stage = FALSE,
                   norm.option = c("all", "by.group"),
                   rescale = c("residSD", "none", "model"),
                   start = NULL,
                   AIREML.tol = 1e-4,
                   max.iter = 100,
                   EM.iter = 0,
                   drop.zeros = TRUE,
                   return.small = FALSE,
                   verbose = TRUE) {
              
              family <- .checkFamily(family, two.stage)
              
              if (is.data.table(x)) x <- as.data.frame(x)

              desmat <- createDesignMatrix(x, outcome, covars, group.var)
              # if there was missing data, need to subset cov.mat
              if (!is.null(cov.mat)) {
                  .checkRownames(cov.mat, x)
                  if (nrow(desmat$X) < nrow(x)) {
                      ind <- which(rownames(x) %in% rownames(desmat$X))
                      cov.mat <- .covMatSubset(cov.mat, ind)
                  }
                  cov.mat <-  .setCovMatNames(cov.mat)
              }


              null.model <- .fitNullModel(y=desmat$y, X=desmat$X, covMatList=cov.mat,
                                          group.idx=desmat$group.idx, family=family,
                                          start=start, AIREML.tol=AIREML.tol,
                                          max.iter=max.iter, EM.iter=EM.iter,
                                          drop.zeros=drop.zeros,
                                          return.small=return.small, verbose=verbose)
              rownames(null.model$fit) <- rownames(desmat$y)

              # Add model string elements here because we need the outcome string.
              null.model$model$outcome <- .modelOutcomeString(outcome, inverse_normal = FALSE)
              null.model$model$covars <- covars
              null.model$model$formula <- .modelString(outcome, covars = covars, random = names(cov.mat),
                                                       group.var = group.var, inverse_normal = FALSE)

              if(two.stage){
                  # fit the second stage model
                  null.model <- nullModelInvNorm(null.model, cov.mat=cov.mat, 
                                                 norm.option=norm.option, rescale=rescale, 
                                                 AIREML.tol=AIREML.tol, max.iter=max.iter,
                                                 EM.iter=EM.iter, drop.zeros=drop.zeros, 
                                                 return.small=return.small, verbose=verbose)
              }

              null.model
          })

setMethod("fitNullModel",
          "AnnotatedDataFrame",
          function(x, outcome,
                   covars = NULL,
                   cov.mat = NULL,
                   group.var = NULL,
                   sample.id = NULL,
                   ...) {

              x <- pData(x)
              if (is(x, "tbl")) x <- as.data.frame(x)
              rownames(x) <- x$sample.id

              if (!is.null(cov.mat)) {
                  .checkSampleId(cov.mat, x)
              }

              ## subset data.frame and cov.mat for selected samples
              if (!is.null(sample.id)) {
                  stopifnot(all(sample.id %in% x$sample.id))
                  ind <- x$sample.id %in% sample.id
                  x <- x[ind,]
                  if (!is.null(cov.mat)) {
                      ind <- which(.covMatNames(cov.mat) %in% sample.id)
                      cov.mat <- .covMatSubset(cov.mat, ind)
                  }
              }

              ## reorder data.frame to match cov.mat
              if (!is.null(cov.mat)) {
                  ind <- match(.covMatNames(cov.mat), rownames(x))
                  x <- x[ind,]
              }

              nm <- fitNullModel(x, outcome, covars, cov.mat, group.var, ...)
              # XXX: Decide if we should reorder columns in the fit data frame such that sample.id is the first column.
              nm$fit$sample.id <- rownames(nm$model.matrix)
              nm
          })

setMethod("fitNullModel",
          "SeqVarData",
          function(x, ...) {
              fitNullModel(sampleData(x), ...)
          })

setMethod("fitNullModel",
          "ScanAnnotationDataFrame",
          function(x, ...) {
              class(x) <- "AnnotatedDataFrame"
              varLabels(x)[varLabels(x) == "scanID"] <- "sample.id"
              fitNullModel(x, ...)
          })

setMethod("fitNullModel",
          "GenotypeData",
          function(x, ...) {
              fitNullModel(getScanAnnotation(x), ...)
          })


## function to check 'family' input
.checkFamily <- function(family, two.stage=FALSE) {
    
    if(is.character(family)){
        family <- get(family)
    }
    if(is.function(family)){
        family <- family()
    }
    if(is.null(family$family)){
        stop("'family' not recognized")
    }
    if (!is.element(family$family, c("gaussian", "binomial", "poisson"))){
        stop("family must be one of gaussian, binomial, or poisson")
    }
    
    if((two.stage) & (family$family != "gaussian")){
        stop('two stage model only applies when family is "gaussian"')
    }
    
    return(family)
}


## function to get sample.ids from cov.mat rownames and/or column names
## error if they don't match
.covMatNames <- function(cov.mat) {
    if (is.list(cov.mat) & length(cov.mat) == 1) {
        cov.mat <- cov.mat[[1]]
    }
    if (!is.list(cov.mat)) {
        names <- dimnames(cov.mat)
        if (!is.null(names[[1]]) & !is.null(names[[2]]) & any(names[[1]] != names[[2]])) {
            stop("dimnames of cov.mat should be identical")
        }
        if (!is.null(names[[1]])) {
            return(names[[1]])
        } else {
            return(names[[2]])
        }
    } else {
        rows <- lapply(cov.mat, rownames)
        cols <- lapply(cov.mat, colnames)
        if (!.listIdentical(rows)) stop("dimnames of cov.mat should be identical")
        if (!.listIdentical(cols)) stop("dimnames of cov.mat should be identical")
        if (!is.null(rows[[1]])) {
            return(rows[[1]])
        } else {
            return(cols[[1]])
        }
    }
}


## need to subset cov.mat in case of missing data
## don't re-order in case we have a block diagonal matrix
.covMatSubset <- function(cov.mat, index) {
    .subset <- function(x, index) {
        if (identical(index, 1:nrow(x))) {
            return(x)
        } else {
            return(x[index,index])
        }
    }
    if (!is.list(cov.mat)) {
        return(.subset(cov.mat, index))
    } else {
        return(lapply(cov.mat, .subset, index))
    }
}


## match rownames between cov.mat and data frame
.checkRownames <- function(cov.mat, x) {
    nms <- .covMatNames(cov.mat)
    if (!is.null(nms) & !all(nms == rownames(x))) {
        stop("dimnames of cov.mat must match rownames of x")
    }
}


## match sample.id between cov.mat and data frame
.checkSampleId <- function(cov.mat, x) {
    nms <- .covMatNames(cov.mat)
    if (is.null(nms)) {
        stop("provide sample.id in rownames and/or colnames of cov.mat")
    } else if (!all(nms %in% x$sample.id)) {
        stop("all sample names in dimnames of cov.mat must be present in x$sample.id")
    }
}


## Return a small version of the null model (no NxN matrices)
nullModelSmall <- function(null.model) {
    null.model$cholSigmaInv <- NULL
    null.model$CX <- NULL
    null.model$CXCXI <- NULL
    null.model
}


isNullModelSmall <- function(null.model) {
    is.null(null.model$cholSigmaInv)
}


.modelString <- function(outcome, covars = NULL, random = NULL, group.var = NULL, inverse_normal = FALSE) {
    model.outcome <- .modelOutcomeString(outcome, inverse_normal = inverse_normal)
    model.covars <- if (is.null(covars)) NULL else .modelCovarString(covars)
    model.random <- if (is.null(random)) NULL else paste(paste0("(1|", random, ")"), collapse=" + ")
    model.var <- if (is.null(group.var)) NULL else paste0("var(", group.var, ")")
    model.string <- paste(c(model.covars, model.random, model.var), collapse=" + ")
    paste(model.outcome, model.string, sep=" ~ ")
}

.modelCovarString <- function(covars) {
  paste(covars, collapse=" + ")
}

.modelOutcomeString <- function(outcome, inverse_normal = FALSE) {
  if (inverse_normal) {
    return(sprintf("rankInvNorm(resid(%s))", outcome))
  } else {
    return(outcome)
  }
}

.setCovMatNames <- function(cov.mat) {
  if (!is.list(cov.mat)) {
    cov.mat <- list(cov.mat)
  }
  # What about the case where one element of the matrix has names, and the others don't?
  if (is.null(names(cov.mat))) {
    names(cov.mat) <- LETTERS[1:length(cov.mat)]
  } else {
    if (any(names(cov.mat) == "")) {
      stop("Some names for cov.mat list are missing.")
    }
  }
  # Check that no names are duplicated.
  if (any(duplicated(names(cov.mat)))) {
    stop("Some names for cov.mat list are duplicated.")
  }
  cov.mat
}

# Update a null model from the old format to the new format.
.updateNullModelFormat <- function(nullmod) {

  msg <- paste(
    "This null model was created with an older version of GENESIS and is being updated to use the current version.",
    "Model formula strings may not exist.",
    "Please consider re-running fitNullModel with the current GENESIS version."
  )

  # Check if old format or new format.
  if (!("fit" %in% names(nullmod))) {

    warning(msg)

    outcome_name <- colnames(nullmod$outcome)

    # Update.
    nullmod$fit <- data.frame(
      outcome = as.vector(nullmod$outcome),
      workingY = as.vector(nullmod$workingY),
      fitted.values = as.vector(nullmod$fitted.values),
      resid.marginal = as.vector(nullmod$resid.marginal),
      resid.PY = as.vector(nullmod$resid),
      resid.cholesky = as.vector(nullmod$Ytilde),
      stringsAsFactors = FALSE
    )
    optional_elements <- c("resid.conditional", "sample.id")
    for (element in optional_elements) {
      if (element %in% names(nullmod)) {
        nullmod$fit[[element]] <- nullmod[[element]]
        nullmod[[element]] <- NULL
      }
    }
    rownames(nullmod$fit) <- rownames(nullmod$model.matrix)

    # Add linear.predictor for mixed models.
    if (nullmod$family$mixedmodel) {
      nullmod$fit$linear.predictor <- nullmod$fit$workingY - nullmod$fit$resid.conditional
    } else {
      # Remove resid.conditional if it exists and the model is not a mixed model.
      # This happened in the past because resid.conditional was calculated for WLS models, but didn't need to be.
      nullmod$fit$resid.conditional <- NULL
    }

    # Remove elements that were moved into fit.
    nullmod$outcome <- NULL
    nullmod$workingY <- NULL
    nullmod$fitted.values <- NULL
    nullmod$resid.marginal <- NULL
    nullmod$resid <- NULL
    nullmod$Ytilde <- NULL

    # Add model element. This element was added along with the fit data frame,
    # so there shouldn't be a case where there is a "fit" element but there is
    # not a "model" element.
    nullmod$model <- list(
      family = nullmod$family,
      hetResid = nullmod$hetResid,
      outcome = outcome_name
    )
    nullmod$family <- NULL
    nullmod$hetResid <- NULL
  }

  # Add RSS0
  if (is.null(nullmod$RSS0)) {
    # This is a minor update - don't need to print the warning for it.
    nullmod$RSS0 <- as.numeric(crossprod(nullmod$fit$resid.cholesky))
  }

  # Add W matrix for mixed models
  if(nullmod$model$family$mixedmodel){
    if(nullmod$model$family$family == "gaussian"){
      n <- nrow(nullmod$fit)
      m <- length(nullmod$varComp) - length(nullmod$group.idx)
      if(is.null(nullmod$group.idx)){
        diagV <- rep(nullmod$varComp[m+1], n)
      }else{
        g <- length(nullmod$group.idx)
        diagV <- rep(NA, n)
        for(i in 1:g){
          diagV[nullmod$group.idx[[i]]] <- nullmod$varComp[m+i]
        }
      }
      nullmod$W <- 1/diagV
    }else{
      eta <- nullmod$fit$linear.predictor
      mu <- nullmod$model$family$linkinv(eta)
      vmu <- nullmod$model$family$variance(mu)
      gmuinv <- nullmod$model$family$mu.eta(eta)
      diagV <- as.vector(vmu)/as.vector(gmuinv)^2
      nullmod$W <- 1/diagV
    }
  }

  nullmod
}
