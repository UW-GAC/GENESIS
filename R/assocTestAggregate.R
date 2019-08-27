setGeneric("assocTestAggregate", function(gdsobj, ...) standardGeneric("assocTestAggregate"))

.match.arg <- function(test) {
    if (length(test) > 1) test <- NULL
    match.arg(test, choices=c("Burden", "SKAT", "fastSKAT", "SMMAT", "fastSMMAT", "SKATO"))
}

setMethod("assocTestAggregate",
          "SeqVarIterator",
          function(gdsobj, null.model, AF.max=1,
                   weight.beta=c(1,1), weight.user=NULL,
                   test=c("Burden", "SKAT", "fastSKAT", "SMMAT", "SKATO"),
                   burden.test=c("Score", "Wald"),
                   # pval.method=c("davies", "kuonen", "liu"),
                   neig = 200, ntrace = 500,
                   rho = seq(from = 0, to = 1, by = 0.1),
                   sparse=TRUE, imputed=FALSE,
                   male.diploid=TRUE, genome.build=c("hg19", "hg38"),
                   verbose=TRUE) {

              # check argument values
              test <- .match.arg(test)
              burden.test <- match.arg(burden.test)
              # pval.method <- match.arg(pval.method)

              # don't use sparse matrices for imputed dosages
              if (imputed) sparse <- FALSE

              # coerce null.model if necessary
              if (sparse) null.model <- .nullModelAsMatrix(null.model)
              
              # filter samples to match null model
              sample.index <- .setFilterNullModel(gdsobj, null.model, verbose=verbose)

              # do we need to match on alleles?
              match.alleles <- any(c("ref", "alt") %in% names(mcols(currentRanges(gdsobj))))

              # check ploidy
              if (SeqVarTools:::.ploidy(gdsobj) == 1) male.diploid <- FALSE
              
              # results
              res <- list()
              res.var <- list()
              i <- 1
              n.iter <- length(variantFilter(gdsobj))
              set.messages <- ceiling(n.iter / 100) # max messages = 100
              iterate <- TRUE
              while (iterate) {
                  var.info <- variantInfo(gdsobj, alleles=match.alleles, expanded=TRUE)
                  
                  if (!imputed) {
                      geno <- expandedAltDosage(gdsobj, use.names=FALSE, sparse=sparse)[sample.index,,drop=FALSE]
                  } else {
                      geno <- imputedDosage(gdsobj, use.names=FALSE)[sample.index,,drop=FALSE]
                  }

                  if (match.alleles) {
                      index <- .matchAlleles(gdsobj, var.info)
                      var.info <- var.info[index,,drop=FALSE]
                      geno <- geno[,index,drop=FALSE]
                  } else {
                      index <- NULL
                  }

                  # number of non-missing samples
                  # n.obs <- colSums(!is.na(geno))
                  n.obs <- .countNonMissing(geno, MARGIN = 2)
                  
                  # allele frequency
                  freq <- .alleleFreq(gdsobj, geno, variant.index=index, sample.index=sample.index,
                                      male.diploid=male.diploid, genome.build=genome.build)
                  
                  # filter monomorphic variants
                  keep <- .filterMonomorphic(geno, count=n.obs, freq=freq$freq, imputed=imputed)

                  # exclude variants with freq > max
                  keep <-  keep & freq$freq <= AF.max
                  if (!all(keep)) {
                      var.info <- var.info[keep,,drop=FALSE]
                      geno <- geno[,keep,drop=FALSE]
                      n.obs <- n.obs[keep]
                      freq <- freq[keep,,drop=FALSE]
                  }

                  # weights
                  if (is.null(weight.user)) {
                      # Beta weights
                      weight <- .weightFromFreq(freq$freq, weight.beta)
                  } else {
                      # user supplied weights
                      weight <- currentVariants(gdsobj)[[weight.user]][expandedVariantIndex(gdsobj)]
                      if (!is.null(index)) weight <- weight[index]
                      weight <- weight[keep]
                      
                      weight0 <- is.na(weight) | weight == 0
                      if (any(weight0)) {
                          keep <- !weight0
                          var.info <- var.info[keep,,drop=FALSE]
                          geno <- geno[,keep,drop=FALSE]
                          n.obs <- n.obs[keep]
                          freq <- freq[keep,,drop=FALSE]
                          weight <- weight[keep]
                      }
                  }
                  
                  # number of variant sites
                  n.site <- length(unique(var.info$variant.id))

                  # number of alternate alleles
                  n.alt <- sum(geno, na.rm=TRUE)
                  
                  # number of samples with observed alternate alleles > 0
                  n.sample.alt <- sum(rowSums(geno, na.rm=TRUE) > 0)
               
                  res[[i]] <- data.frame(n.site, n.alt, n.sample.alt)
                  res.var[[i]] <- cbind(var.info, n.obs, freq, weight)
                  
                  if (n.site > 0) {
                      # mean impute missing values
                      if (any(n.obs < nrow(geno))) {
                          geno <- .meanImpute(geno, freq$freq)
                      }

                      # do the test
                      assoc <- testVariantSet(null.model, G=geno, weights=weight, 
                                              test=test, burden.test=burden.test, 
                                              neig = neig, ntrace = ntrace,
                                              rho=rho)
                                              # pval.method=pval.method)
                      res[[i]] <- cbind(res[[i]], assoc, stringsAsFactors=FALSE)
                  }

                  if (verbose & n.iter > 1 & i %% set.messages == 0) {
                      message(paste("Iteration", i , "of", n.iter, "completed"))
                  }
                  i <- i + 1
                  iterate <- iterateFilter(gdsobj, verbose=FALSE)
              }

              res <- list(results=bind_rows(res), variantInfo=res.var)
              .annotateAssoc(gdsobj, res)
          })



setMethod("assocTestAggregate",
          "GenotypeIterator",
          function(gdsobj, null.model, AF.max=1,
                   weight.beta=c(1,1), weight.user=NULL,
                   test=c("Burden", "SKAT", "fastSKAT", "SMMAT", "SKATO"),
                   burden.test=c("Score", "Wald"),
                   # pval.method=c("davies", "kuonen", "liu"),
                   neig = 200, ntrace = 500,
                   rho = seq(from = 0, to = 1, by = 0.1),
                   male.diploid=TRUE, verbose=TRUE) {

              # check argument values
              test <- .match.arg(test)
              burden.test <- match.arg(burden.test)
              # pval.method <- match.arg(pval.method)

              # filter samples to match null model
              sample.index <- .sampleIndexNullModel(gdsobj, null.model)

              # results
              res <- list()
              res.var <- list()
              i <- 1
              n.iter <- length(snpFilter(gdsobj))
              set.messages <- ceiling(n.iter / 100) # max messages = 100
              iterate <- TRUE
              while (iterate) {
                  var.info <- variantInfo(gdsobj)
                  
                  geno <- getGenotypeSelection(gdsobj, scan=sample.index, order="selection",
                                               transpose=TRUE, use.names=FALSE, drop=FALSE)
                  
                  # number of non-missing samples
                  # n.obs <- colSums(!is.na(geno))
                  n.obs <- .countNonMissing(geno, MARGIN = 2)
                  
                  # allele frequency
                  freq <- .alleleFreq(gdsobj, geno, sample.index=sample.index,
                                      male.diploid=male.diploid)
                  
                  # filter monomorphic variants
                  keep <- .filterMonomorphic(geno, count=n.obs, freq=freq$freq)

                  # exclude variants with freq > max
                  keep <-  keep & freq$freq <= AF.max
                  if (!all(keep)) {
                      var.info <- var.info[keep,,drop=FALSE]
                      geno <- geno[,keep,drop=FALSE]
                      n.obs <- n.obs[keep]
                      freq <- freq[keep,,drop=FALSE]
                  }

                  # weights
                  if (is.null(weight.user)) {
                      # Beta weights
                      weight <- .weightFromFreq(freq$freq, weight.beta)
                  } else {
                      # user supplied weights
                      weight <- getSnpVariable(gdsobj, weight.user)
                      weight <- weight[keep]
                      
                      weight0 <- is.na(weight) | weight == 0
                      if (any(weight0)) {
                          keep <- !weight0
                          var.info <- var.info[keep,,drop=FALSE]
                          geno <- geno[,keep,drop=FALSE]
                          n.obs <- n.obs[keep]
                          freq <- freq[keep,,drop=FALSE]
                          weight <- weight[keep]
                      }
                  }
                  
                  # number of variant sites
                  n.site <- length(unique(var.info$variant.id))

                  # number of alternate alleles
                  n.alt <- sum(geno, na.rm=TRUE)
                  
                  # number of samples with observed alternate alleles > 0
                  n.sample.alt <- sum(rowSums(geno, na.rm=TRUE) > 0)
               
                  res[[i]] <- data.frame(n.site, n.alt, n.sample.alt)
                  res.var[[i]] <- cbind(var.info, n.obs, freq, weight)
                  
                  if (n.site > 0) {
                      # mean impute missing values
                      if (any(n.obs < nrow(geno))) {
                          geno <- .meanImpute(geno, freq$freq)
                      }

                      # do the test
                      assoc <- testVariantSet(null.model, G=geno, weights=weight, 
                                              test=test, burden.test=burden.test, 
                                              neig = neig, ntrace = ntrace,
                                              rho=rho)
                                              # pval.method=pval.method)
                      res[[i]] <- cbind(res[[i]], assoc, stringsAsFactors=FALSE)
                  }

                  if (verbose & n.iter > 1 & i %% set.messages == 0) {
                      message(paste("Iteration", i , "of", n.iter, "completed"))
                  }
                  i <- i + 1
                  iterate <- GWASTools::iterateFilter(gdsobj)
              }

              res <- list(results=bind_rows(res), variantInfo=res.var)
              .annotateAssoc(gdsobj, res)
          })
