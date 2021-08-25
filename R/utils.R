# GenotypeData methods for SeqVarTools generics
setMethod("chromWithPAR",
          "GenotypeData",
          function(gdsobj, ...) {
              getChromosome(gdsobj, char=TRUE)
          })

setMethod("validateSex",
          "GenotypeData",
          function(x) {
              sex <- suppressWarnings(getSex(x))
              if (!is.null(sex)) {
                  if (all(sex %in% c(1,2,NA))) {
                      sex <- c("M", "F")[as.integer(sex)]
                  }
                  if (!all(sex %in% c("M", "F", NA))) {
                      sex <- NULL
                  }
              }
              sex
          })

setMethod("variantInfo",
          "GenotypeData",
          function(gdsobj, alleles=FALSE) {
              data.frame(variant.id=getSnpID(gdsobj),
                         chr=getChromosome(gdsobj, char=TRUE),
                         pos=getPosition(gdsobj),
                         stringsAsFactors=FALSE)
          })

setMethod("variantFilter",
          "GenotypeIterator",
          function(x) {
              snpFilter(x)
          })


# function to pre-process genotype data before testing
.prepGenoBlock <- function(x, AF.max=1, geno.coding="additive", imputed=FALSE, 
                           sex=NULL, male.diploid=TRUE) {
    
    var.info <- x$var.info
    geno <- x$geno
    chr <- x$chr
    weight <- x$weight # only applies to aggregate tests, NULL otherwise
    rm(x)
    
    # take note of number of non-missing samples
    #n.obs <- colSums(!is.na(geno))
    n.obs <- .countNonMissing(geno, MARGIN = 2)
    
    # allele frequency
    freq <- .alleleFreq(geno, chr, sex, male.diploid=male.diploid)
    
    # filter monomorphic variants
    keep <- .filterMonomorphic(geno, count=n.obs, freq=freq$freq, imputed=imputed)
    
    # exclude variants with freq > max
    keep <-  keep & freq$freq <= AF.max
    
    # exclude variants with weight 0
    if (!is.null(weight)) {
        weight0 <- is.na(weight) | weight == 0
        if (any(weight0)) {
            keep <- keep & !weight0
        }
    }
    
    # recessive or dominant coding
    if (geno.coding != "additive") {
        if (geno.coding == "recessive") {
            ## if wanting to test a recessive model, the genotypes 0 and 1 are '0'
            ## and the genotype 2 is '1',
            ## e.g. indicator of participant having 2 copies of the alternate allele
            geno[geno == 1] <- 0L
            geno[geno == 2] <- 1L
            out.col <- "n.hom.alt"
        } else if (geno.coding == "dominant") {
            geno[geno == 2] <- 1L
            out.col <- "n.any.alt"
        }
        
        # count number of carriers
        freq[[out.col]] <- colSums(geno, na.rm=TRUE)
        
        # remove variants which are now monomorphic (all 0s or all 1s)
        keep <- keep & !(freq[[out.col]] == 0 | freq[[out.col]] == n.obs)
    }
    
    if (!all(keep)) {
        var.info <- var.info[keep,,drop=FALSE]
        geno <- geno[,keep,drop=FALSE]
        n.obs <- n.obs[keep]
        freq <- freq[keep,,drop=FALSE]
        weight <- weight[keep]
    }
    
    return(list(var.info=var.info, n.obs=n.obs, freq=freq, geno=geno, weight=weight))
}


# check for all genotypes identical (including all hets)
# allow a tolerance in case we have imputed dosage values rather than integer
# return an index for subsetting rather than modifying geno
# (so we can subset variant info also)
# we have to compute count and freq anyway, so pass as arguments
.filterMonomorphic <- function(geno, count, freq, imputed=FALSE) {
    #count <- colSums(!is.na(geno))
    if (!imputed) {
        #freq <- 0.5*colMeans(geno, na.rm=TRUE)
        isref <- freq == 0
        isalt <- freq == 1
        ishet <- colSums(geno == 1, na.rm=TRUE) == count
    } else {
        tol <- .Machine$double.eps ^ 0.5
        isref <- colSums(abs(geno) < tol, na.rm=TRUE) == count
        isalt <- colSums(abs(geno - 1) < tol, na.rm=TRUE) == count
        ishet <- colSums(abs(geno - 2) < tol, na.rm=TRUE) == count
    }
    !(isref | isalt | ishet)
}


.alleleFreq <- function(geno, chr=NULL, sex=NULL, male.diploid=TRUE) {

    if (is.null(sex) | is.null(chr)) {
        #freq <- 0.5*colMeans(geno, na.rm=TRUE)
        count <- colSums(geno, na.rm=TRUE)
        # nsamp <- colSums(!is.na(geno))
        nsamp <- .countNonMissing(geno, MARGIN = 2)
        freq <- count/(2*nsamp)
        mac <- round(pmin(count, 2*nsamp - count))
        return(data.frame(freq=freq, MAC=mac))
    }

    # check chromosome
    X <- chr %in% "X"
    Y <- chr %in% "Y"
    auto <- !X & !Y

    # allele count vectors
    #freq <- rep(NA, ncol(geno))
    count <- rep(NA, ncol(geno))
    possible <- rep(NA, ncol(geno))

    # autosomes
    if (any(auto)) {
        #freq[auto] <- 0.5*colMeans(geno[, auto, drop=FALSE], na.rm=TRUE)
        count[auto] <- colSums(geno[, auto, drop=FALSE], na.rm=TRUE)
        # possible[auto] <- 2 * colSums(!is.na(geno[, auto, drop=FALSE]))
        possible[auto] <- 2 * .countNonMissing(geno[, auto, drop=FALSE], MARGIN = 2)
    }

    # X chrom
    if (any(X)) {
        female <- sex %in% "F"
        male <- sex %in% "M"
        F.count <- colSums(geno[female, X, drop=FALSE], na.rm=TRUE)
        # F.nsamp <- colSums(!is.na(geno[female, X, drop=FALSE]))
        F.nsamp <- .countNonMissing(geno[female, X, drop=FALSE], MARGIN = 2)
        M.count <- colSums(geno[male, X, drop=FALSE], na.rm=TRUE)
        if (male.diploid) {
            M.count <- 0.5*M.count
        }
        # M.nsamp <- colSums(!is.na(geno[male, X, drop=FALSE]))
        M.nsamp <- .countNonMissing(geno[male, X, drop=FALSE], MARGIN = 2)
        #freq[X] <- (F.count + M.count)/(2*F.nsamp + M.nsamp)
        count[X] <- (F.count + M.count)
        possible[X] <- (2*F.nsamp + M.nsamp)
    }

    # Y chrom
    if (any(Y)) {
        male <- sex %in% "M"
        #freq[Y] <- colMeans(geno[male, Y, drop=FALSE], na.rm=TRUE)
        count[Y] <- colSums(geno[male, Y, drop=FALSE], na.rm=TRUE)
        if (male.diploid) {
            #freq[Y] <- 0.5*freq[Y]
            count[Y] <- 0.5*count[Y]
        }
        # possible[Y] <- colSums(!is.na(geno[male, Y, drop=FALSE]))
        possible[Y] <- .countNonMissing(geno[male, Y, drop=FALSE], MARGIN = 2)
    }

    freq <- count / possible
    mac <- round(pmin(count, possible - count))
    data.frame(freq=freq, MAC=mac)
}


.meanImputeFn <- function(geno, freq, geno.coding="additive") {
    miss.idx <- which(is.na(geno))
    miss.var.idx <- ceiling(miss.idx/nrow(geno))
    if (geno.coding == "additive") {
        geno[miss.idx] <- 2*freq[miss.var.idx]
    } else {
        freq <- colMeans(geno, na.rm=TRUE)
        geno[miss.idx] <- freq[miss.var.idx]
    }
    geno
}


.weightFromFreq <- function(freq, weight.beta) {
    freq <- pmin(freq, 1-freq)
    dbeta(freq, weight.beta[1], weight.beta[2])
}


# set a sample filter, and return the index to put filtered samples
# in the same order as the null model
.setFilterNullModel <- function(gdsobj, null.model, verbose=TRUE) {
    if (!is.null(null.model$fit$sample.id)) {
        seqSetFilter(gdsobj, sample.id=null.model$fit$sample.id, verbose=verbose)
        sample.index <- match(null.model$fit$sample.id, seqGetData(gdsobj, "sample.id"))
    } else {
        sample.index <- seq_along(seqGetData(gdsobj, "sample.id"))
    }
    sample.index
}


# return sample.index to match GenotypeData object to a null model
.sampleIndexNullModel <- function(gdsobj, null.model) {
    if (!is.null(null.model$fit$sample.id)) {
        sample.index <- match(null.model$fit$sample.id, getScanID(gdsobj))
    } else {
        sample.index <- match(rownames(null.model$model.matrix),
                              sampleNames(getScanAnnotation(gdsobj)))
    }
    sample.index
}


.modelMatrixColumns <- function(null.model, col.name) {
    cols <- unlist(lapply(col.name, function(x) grep(paste0("^", x), colnames(null.model$model.matrix))))
    null.model$model.matrix[,cols,drop=FALSE]
}


.matchAlleles <- function(gdsobj, var.info) {
    seqnames <- NULL
    if (nrow(var.info) == 0) return(integer(0))
    var.info$n <- 1:nrow(var.info)
    var.sel <- as.data.table(currentRanges(gdsobj))
    var.sel <- var.sel[,`:=`(seqnames=as.character(seqnames))]
    setnames(var.sel, c("seqnames", "start"), c("chr", "pos"))
    match.cols <- c("chr", "pos")
    if ("ref" %in% names(var.sel)) match.cols <- c(match.cols, "ref")
    if ("alt" %in% names(var.sel)) match.cols <- c(match.cols, "alt")
    var.match <- merge(as.data.table(var.info), var.sel, by=match.cols)
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

.stopOnError <- function(x) {
    err.chk <- sapply(x, is, "error")
    if (any(err.chk)) {
        ind <- which(err.chk)[1]
        # only in serial execution is ind guaranteed to be the right index in an error state
        #message("Error detected in iteration ", ind)
        stop(x[[ind]])
    }
}
