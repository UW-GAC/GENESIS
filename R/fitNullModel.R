
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
              desmat <- createDesignMatrix2(x, outcome, covars, group.var)
              fitNullMod(y=desmat$y, X=desmat$X, covMatList=cov.mat,
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
              desmat <- createDesignMatrix2(x, outcome, covars, group.var, sample.id)

              # subset or re-order cov.mat if necessary
              if (!is.null(cov.mat)) {
                  if (!is.list(cov.mat)) {
                      cov.mat <- list(A=cov.mat)
                  }
                  if (!is.null(sample.id)) {
                      cov.mat <- lapply(cov.mat, .orderSamples,
                                           orig.ids=x$sample.id, new.ids=sample.id)
                  }
              }
              
              nm <- fitNullMod(y=desmat$y, X=desmat$X, covMatList=cov.mat,
                               group.idx=desmat$group.idx, ...)
              nm$sample.id <- rownames(nm$model.matrix)
              nm
          })

setMethod("fitNullModel",
          "SeqVarData",
          function(x, ...) {
              fitNullModel(sampleData(x), ...)
          })


nullModelInvNorm <- function(null.model, cov.mat = NULL,
                             norm.option = c("by.group", "all"),
                             rescale = c("none", "model", "residSD"),
                             AIREML.tol = 1e-6, max.iter = 100, verbose = TRUE) {

    # subset or re-order cov.mat if necessary
    if (!is.null(cov.mat)) {
        if (!is.list(cov.mat)) {
            cov.mat <- list(A=cov.mat)
        }
        if (!is.null(null.model$sample.id)) {
            cov.mat <- lapply(cov.mat, function(x) {
                .orderSamples(x, orig.ids=rownames(x), new.ids=null.model$sample.id)
            })
        }
    }

    updateNullModOutcome(null.model, covMatList=cov.mat, rankNorm.option=norm.option,
                         rescale=rescale, AIREML.tol=AIREML.tol, max.iter=max.iter,
                         verbose=verbose)
}


.orderSamples <- function(cov.mat, orig.ids, new.ids) {
    if (!is.null(rownames(cov.mat)) & !is.null(colnames(cov.mat))) {
        stopifnot(identical(rownames(cov.mat), colnames(cov.mat)))
        stopifnot(all(new.ids %in% rownames(cov.mat)))
    } else if (!is.null(rownames(cov.mat))) {
        stopifnot(all(new.ids %in% rownames(cov.mat)))
        colnames(cov.mat) <- rownames(cov.mat)
    } else if (!is.null(colnames(cov.mat))) {
        stopifnot(all(new.ids %in% colnames(cov.mat)))
        rownames(cov.mat) <- colnames(cov.mat)
    } else {
        warning("no dimnames given for cov.mat; assuming order of samples matches data frame")
        dimnames(cov.mat) <- list(orig.ids, orig.ids)
    }
    orig.ids <- as.character(orig.ids)
    new.ids <- as.character(new.ids)
    if (identical(rownames(cov.mat), new.ids)) {
        return(cov.mat)
    } else {
        return(cov.mat[new.ids, new.ids])
    }
}
