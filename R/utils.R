# GenotypeData methods for SeqVarTools generics
setMethod("chromWithPAR",
          "GenotypeData",
          function(gdsobj) {
              getChromosome(gdsobj, char=TRUE)
          })

setMethod("validateSex",
          "GenotypeData",
          function(x) {
              sex <- suppressWarnings(getSex(x))
              if (!is.null(sex)) {
                  if (all(sex %in% c(1,2,NA))) {
                      sex <- c("M", "F")[sex]
                  }
                  if (!all(sex %in% c("M", "F", NA))) {
                      sex <- NULL
                  }
              }
              sex
          })

# index is in case we had to subset geno so it no longer matches the variant filter
# (in the case of allele matching)
.alleleFreq <- function(gdsobj, geno, variant.index=NULL, sample.index=NULL) {

    # check sex
    sex <- validateSex(gdsobj)
    if (!is.null(sample.index)) sex <- sex[sample.index]
    if (is.null(sex)) {
        return(0.5*colMeans(geno, na.rm=TRUE))
    }

    # check chromosome
    chr <- chromWithPAR(gdsobj)
    if (!is.null(variant.index)) chr <- chr[variant.index]
    X <- chr %in% "X"
    Y <- chr %in% "Y"
    auto <- !X & !Y

    # allele frequency vector
    freq <- rep(NA, ncol(geno))

    # autosomes
    if (any(auto)) {
        freq[auto] <- 0.5*colMeans(geno[, auto, drop=FALSE], na.rm=TRUE)
    }

    # X chrom
    if (any(X)) {
        female <- sex %in% "F"
        male <- sex %in% "M"
        F.count <- colSums(geno[female, X, drop=FALSE], na.rm=TRUE)
        F.nsamp <- colSums(!is.na(geno[female, X, drop=FALSE]))
        M.count <- 0.5*colSums(geno[male, X, drop=FALSE], na.rm=TRUE)
        M.nsamp <- colSums(!is.na(geno[male, X, drop=FALSE]))
        freq[X] <- (F.count + M.count)/(2*F.nsamp + M.nsamp)
    }

    # Y chrom
    if (any(Y)) {
        male <- sex %in% "M"
        freq[Y] <- 0.5*colMeans(geno[male, Y, drop=FALSE], na.rm=TRUE)
    }

    freq
}

.meanImpute <- function(geno, freq) {
    miss.idx <- which(is.na(geno))
    miss.var.idx <- ceiling(miss.idx/nrow(geno))
    geno[miss.idx] <- 2*freq[miss.var.idx]
    geno
}

.weightFromFreq <- function(freq, weight.beta) {
    freq <- pmin(freq, 1-freq)
    dbeta(freq, weight.beta[1], weight.beta[2])
}

# set a sample filter, and return the index to put filtered samples
# in the same order as the null model
.setFilterNullModel <- function(gdsobj, null.model, verbose=TRUE) {
    if (!is.null(null.model$sample.id)) {
        seqSetFilter(gdsobj, sample.id=null.model$sample.id, verbose=verbose)
        match(null.model$sample.id, seqGetData(gdsobj, "sample.id"))
    } else {
        seq_along(seqGetData(gdsobj, "sample.id"))
    }
}

.modelMatrixColumns <- function(null.model, col.name) {
    null.model$model.matrix[,grep(paste0("^", col.name), colnames(null.model$model.matrix)),drop=FALSE]
}

.matchAlleles <- function(gdsobj, var.info) {
    if (nrow(var.info) == 0) return(integer(0))
    var.info$n <- 1:nrow(var.info)
    var.sel <- as.data.frame(currentRanges(gdsobj))
    var.sel$seqnames <- as.character(var.sel$seqnames)
    match.cols <- c("chr"="seqnames", "pos"="start")
    if ("ref" %in% names(var.sel)) match.cols <- c(match.cols, "ref"="ref")
    if ("alt" %in% names(var.sel)) match.cols <- c(match.cols, "alt"="alt")
    var.match <- inner_join(var.info, var.sel, by=match.cols)
    unique(var.match$n)
}


.addNames <- function(iterator, x) {
    gr <- variantRanges(iterator)
    rownames(x$results) <- names(gr)
    names(x$variantInfo) <- names(gr)
    x
}

.addWindows <- function(iterator, x) {
    windows <- variantRanges(iterator)
    win.df <- data.frame(chr=as.character(GenomicRanges::seqnames(windows)),
                         start=BiocGenerics::start(windows),
                         end=BiocGenerics::end(windows),
                         stringsAsFactors=FALSE)
    x$results <- cbind(win.df, x$results)
    x
}

setGeneric(".annotateAssoc", function(gdsobj, x) standardGeneric(".annotateAssoc"))

setMethod(".annotateAssoc",
          "SeqVarWindowIterator",
          function(gdsobj, x) {
              .addWindows(gdsobj, x)
          })

setMethod(".annotateAssoc",
          "SeqVarRangeIterator",
          function(gdsobj, x) {
              .addNames(gdsobj, x)
          })

setMethod(".annotateAssoc",
          "SeqVarListIterator",
          function(gdsobj, x) {
              .addNames(gdsobj, x)
          })

setMethod(".annotateAssoc",
          "SeqVarBlockIterator",
          function(gdsobj, x) {
              x
          })

setMethod(".annotateAssoc",
          "GenotypeIterator",
          function(gdsobj, x) {
              filt.names <- names(snpFilter(gdsobj))
              if (!is.null(filt.names)) {
                  rownames(x$results) <- filt.names
                  names(x$variantInfo) <- filt.names
              }
              x
          })

.listIdentical <- function(x) {
    all(sapply(x[-1], FUN = identical, x[[1]]))
}
