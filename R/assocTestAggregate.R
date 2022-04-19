setGeneric("assocTestAggregate", function(gdsobj, ...) standardGeneric("assocTestAggregate"))

.match.arg <- function(test) {
    if (length(test) > 1) test <- NULL
    match.arg(test, choices=c("Burden", "SKAT", "fastSKAT", "SMMAT", "fastSMMAT", "SKATO", "BinomiRare", "CMP"))
}

setMethod("assocTestAggregate",
          "SeqVarIterator",
          function(gdsobj, null.model, AF.max=1,
                   weight.beta=c(1,1), weight.user=NULL,
                   test=c("Burden", "SKAT", "fastSKAT", "SMMAT", "fastSMMAT", "SKATO", "BinomiRare", "CMP"),
                   neig = 200, ntrace = 500,
                   rho = seq(from = 0, to = 1, by = 0.1),
                   sparse=TRUE, imputed=FALSE,
                   male.diploid=TRUE, genome.build=c("hg19", "hg38"),
                   BPPARAM=bpparam(), verbose=TRUE) {

              # check argument values
              test <- .match.arg(test)

              # for BinomiRare and CMP, restrict to variants where the alternate allele is minor
              if (test %in% c("BinomiRare", "CMP") && AF.max > 0.5) {
                  AF.max <- 0.5
              }

              # don't use sparse matrices for imputed dosages
              if (imputed) sparse <- FALSE

              # Convert old null model format if necessary.
              null.model <- .updateNullModelFormat(null.model)

              # coerce null.model if necessary
              if (sparse) null.model <- .nullModelAsMatrix(null.model)

              # filter samples to match null model
              sample.index <- .setFilterNullModel(gdsobj, null.model, verbose=verbose)

              # get sex for calculating allele freq
              sex <- validateSex(gdsobj)[sample.index]

              # do we need to match on alleles?
              match.alleles <- any(c("ref", "alt") %in% names(mcols(currentRanges(gdsobj))))

              # check ploidy
              if (SeqVarTools:::.ploidy(gdsobj) == 1) male.diploid <- FALSE

              # results
              # n.iter <- length(variantFilter(gdsobj))
              # set.messages <- ceiling(n.iter / 100) # max messages = 100

              if(verbose) message('Using ', bpnworkers(BPPARAM), ' CPU cores')

              i <- 1
              ITER <- function() {
                  iterate <- TRUE
                  if (i > 1) {
                      iterate <- iterateFilter(gdsobj, verbose=FALSE)
                  }
                  i <<- i + 1
                  if (!iterate) {
                      return(NULL)
                  }

                  var.info <- variantInfo(gdsobj, alleles=match.alleles, expanded=TRUE)

                  if (!imputed) {
                      geno <- expandedAltDosage(gdsobj, use.names=FALSE, sparse=sparse)[sample.index,,drop=FALSE]
                  } else {
                      geno <- imputedDosage(gdsobj, use.names=FALSE)[sample.index,,drop=FALSE]
                  }

                  chr <- chromWithPAR(gdsobj, genome.build=genome.build)

                  if (!is.null(weight.user)) {
                      weight <- currentVariants(gdsobj)[[weight.user]][expandedVariantIndex(gdsobj)]
                  } else {
                      weight <- NULL
                  }

                  if (match.alleles) {
                      index <- .matchAlleles(gdsobj, var.info)
                      var.info <- var.info[index,,drop=FALSE]
                      geno <- geno[,index,drop=FALSE]
                      chr <- chr[index]
                      if (!is.null(weight)) {
                          weight <- weight[index]
                      }
                  }


                  return(list(var.info=var.info, geno=geno, chr=chr, weight=weight))
              }

              res <- bpiterate(ITER, .testGenoBlockAggregate, BPPARAM=BPPARAM,
                               sex=sex, null.model=null.model, AF.max=AF.max,
                               weight.beta=weight.beta, test=test,
                               neig=neig, ntrace=ntrace, rho=rho,
                               sparse=sparse, imputed=imputed,
                               male.diploid=male.diploid)
              .stopOnError(res)

              res1 <- lapply(res, function(x) x[[1]])
              res.var <- lapply(res, function(x) x[[2]])
              rm(res)
              res <- list(results=as.data.frame(rbindlist(res1, use.names=TRUE, fill=TRUE)),
                          variantInfo=res.var)
              .annotateAssoc(gdsobj, res)
          })



setMethod("assocTestAggregate",
          "GenotypeIterator",
          function(gdsobj, null.model, AF.max=1,
                   weight.beta=c(1,1), weight.user=NULL,
                   test=c("Burden", "SKAT", "fastSKAT", "SMMAT", "fastSMMAT", "SKATO", "BinomiRare", "CMP"),
                   neig = 200, ntrace = 500,
                   rho = seq(from = 0, to = 1, by = 0.1),
                   male.diploid=TRUE,
                   BPPARAM=bpparam(), verbose=TRUE) {

              # check argument values
              test <- .match.arg(test)

              # for BinomiRare and CMP, restrict to variants where the alternate allele is minor
              if (test %in% c("BinomiRare", "CMP") && AF.max > 0.5) {
                  AF.max  <-  0.5
              }

              # Convert old null model format if necessary.
              null.model <- .updateNullModelFormat(null.model)

              # filter samples to match null model
              sample.index <- .sampleIndexNullModel(gdsobj, null.model)

              # get sex for calculating allele freq
              sex <- validateSex(gdsobj)[sample.index]

              # results
              # n.iter <- length(variantFilter(gdsobj))
              # set.messages <- ceiling(n.iter / 100) # max messages = 100

              if(verbose) message('Using ', bpnworkers(BPPARAM), ' CPU cores')

              i <- 1
              ITER <- function() {
                  iterate <- TRUE
                  if (i > 1) {
                      iterate <- GWASTools::iterateFilter(gdsobj)
                  }
                  i <<- i + 1
                  if (!iterate) {
                      return(NULL)
                  }

                  var.info <- variantInfo(gdsobj)

                  geno <- getGenotypeSelection(gdsobj, scan=sample.index, order="selection",
                                               transpose=TRUE, use.names=FALSE, drop=FALSE)

                  if (!is.null(weight.user)) {
                      weight <- getSnpVariable(gdsobj, weight.user)
                  } else {
                      weight <- NULL
                  }

                  return(list(var.info=var.info, geno=geno, chr=var.info$chr, weight=weight))
              }


              res <- bpiterate(ITER, .testGenoBlockAggregate, BPPARAM=BPPARAM,
                               sex=sex, null.model=null.model, AF.max=AF.max,
                               weight.beta=weight.beta, test=test,
                               neig=neig, ntrace=ntrace, rho=rho,
                               sparse=FALSE, imputed=FALSE,
                               male.diploid=male.diploid)
              .stopOnError(res)

              res1 <- lapply(res, function(x) x[[1]])
              res.var <- lapply(res, function(x) x[[2]])
              rm(res)
              res <- list(results=as.data.frame(rbindlist(res1, use.names=TRUE, fill=TRUE)),
                          variantInfo=res.var)
              .annotateAssoc(gdsobj, res)
          })


.testGenoBlockAggregate <- function(x, sex, null.model, AF.max,
                                    weight.beta, test,
                                    neig, ntrace, rho,
                                    sparse, imputed, male.diploid, ...) {

    x <- .prepGenoBlock(x, AF.max=AF.max, sex=sex, imputed=imputed, male.diploid=male.diploid)
    var.info <- x$var.info
    n.obs <- x$n.obs
    freq <- x$freq
    geno <- x$geno
    weight <- x$weight
    rm(x)

    # weights
    if (is.null(weight)) {
        # Beta weights
        if (test %in% c("BinomiRare", "CMP")){
            weight <- seq(from=1, to=1, length.out=length(freq$freq))
        } else {
            weight <- .weightFromFreq(freq$freq, weight.beta)
        }
    }

    # number of variant sites
    n.site <- length(unique(var.info$variant.id))

    # number of alternate alleles
    n.alt <- sum(geno, na.rm=TRUE)

    # number of samples with observed alternate alleles > 0
    n.sample.alt <- sum(rowSums(geno, na.rm=TRUE) > 0)

    res.i <- data.frame(n.site, n.alt, n.sample.alt)
    res.var.i <- cbind(var.info, n.obs, freq, weight)

    if (n.site > 0) {
        # mean impute missing values
        if (any(n.obs < nrow(geno))) {
            geno <- .meanImpute(geno, freq$freq)
        }

        # do the test
        assoc <- testVariantSet(null.model, G=geno, weights=weight,
                                test=test,
                                neig = neig, ntrace = ntrace,
                                rho=rho)
        res.i <- cbind(res.i, assoc, stringsAsFactors=FALSE)
    }

    # if (verbose & n.iter > 1 & i %% set.messages == 0) {
    #     message(paste("Iteration", i , "of", n.iter, "completed"))
    # }
    return(list(res.i, res.var.i))
}
