
setGeneric("fitNullModel", function(x, ...) standardGeneric("fitNullModel"))

setMethod("fitNullModel",
          "data.frame",
          function(x, outcome,
                   covars = NULL,
                   cov.mat = NULL,
                   group.var = NULL,
                   family = "gaussian",
                   start = NULL,
                   AIREML.tol = 1e-6,
                   max.iter = 100,
                   drop.zeros = TRUE,
                   verbose = TRUE) {
              
              desmat <- createDesignMatrix(x, outcome, covars, group.var)

              # if there was missing data, need to subset cov.mat
              if (!is.null(cov.mat)) {
                  .checkRownames(cov.mat, x)
                  if (nrow(desmat$X) < nrow(x)) {
                      ind <- which(rownames(x) %in% rownames(desmat$X))
                      cov.mat <- .covMatSubset(cov.mat, ind)
                  }
              }
              
              .fitNullModel(y=desmat$y, X=desmat$X, covMatList=cov.mat,
                            group.idx=desmat$group.idx, family=family,
                            start=start, AIREML.tol=AIREML.tol,
                            max.iter=max.iter, drop.zeros=drop.zeros,
                            verbose=verbose)
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
              nm$sample.id <- rownames(nm$model.matrix)
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


nullModelInvNorm <- function(null.model, cov.mat = NULL,
                             norm.option = c("by.group", "all"),
                             rescale = c("none", "model", "residSD"),
                             AIREML.tol = 1e-6, max.iter = 100, verbose = TRUE) {

    updateNullModOutcome(null.model, covMatList=cov.mat, rankNorm.option=norm.option,
                         rescale=rescale, AIREML.tol=AIREML.tol, max.iter=max.iter,
                         verbose=verbose)
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
        .subset(cov.mat, index)
    } else {
        cov.mat <- return(lapply(cov.mat, .subset, index))
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
