setGeneric("assocTestAggregate", function(gdsobj, ...) standardGeneric("assocTestAggregate"))

setMethod("assocTestAggregate",
          "SeqVarIterator",
          function(gdsobj, null.model, AF.max=1,
                   weight.beta=c(1,1), weight.user=NULL,
                   test=c("Burden", "SKAT", "SMMAT"),
                   burden.test=c("Score", "Wald"), rho=0,
                   pval.method=c("davies", "kuonen", "liu"),
                   verbose=TRUE) {

              # check argument values
              test <- match.arg(test)
              burden.test <- match.arg(burden.test)
              pval.method <- match.arg(pval.method)

              # filter samples to match null model
              sample.index <- .setFilterNullModel(gdsobj, null.model, verbose=verbose)

              # do we need to match on alleles?
              match.alleles <- any(c("ref", "alt") %in% names(mcols(currentRanges(gdsobj))))

              # results
              res <- list()
              res.var <- list()
              i <- 1
              n.iter <- length(variantFilter(gdsobj))
              iterate <- TRUE
              while (iterate) {
                  var.info <- variantInfo(gdsobj, alleles=match.alleles, expanded=TRUE)
                  
                  geno <- expandedAltDosage(gdsobj, use.names=FALSE, sparse=TRUE)[sample.index,,drop=FALSE]

                  if (match.alleles) {
                      index <- .matchAlleles(gdsobj, var.info)
                      var.info <- var.info[index,,drop=FALSE]
                      geno <- geno[,index,drop=FALSE]
                  } else {
                      index <- NULL
                  }

                  # allele frequency
                  freq <- .alleleFreq(gdsobj, geno, variant.index=index, sample.index=sample.index)
                  # exclude monomorphic variants
                  mono <- freq %in% c(0,1)
                  # exclude variants with freq > max
                  excl <-  mono | freq > AF.max | is.na(freq)
                  if (any(excl)) {
                      var.info <- var.info[!excl,,drop=FALSE]
                      geno <- geno[,!excl,drop=FALSE]
                      freq <- freq[!excl]
                  }

                  # number of variant sites
                  n.site <- length(unique(var.info$variant.id))

                  # number of alternate alleles
                  n.alt <- sum(geno, na.rm=TRUE)
                  
                  # number of samples with observed alternate alleles > 0
                  n.sample.alt <- sum(rowSums(geno, na.rm=TRUE) > 0)
               
                  # number of non-missing samples
                  n.obs <- colSums(!is.na(geno))
                  
                  # weights
                  if (is.null(weight.user)) {
                      # Beta weights
                      weight <- .weightFromFreq(freq, weight.beta)
                  } else {
                      # user supplied weights
                      weight <- currentVariants(gdsobj)[[weight.user]][expandedVariantIndex(gdsobj)]
                      if (!is.null(index)) weight <- weight[index]
                      weight <- weight[!excl]
                  }
                  
                  res[[i]] <- data.frame(n.site, n.alt, n.sample.alt)
                  res.var[[i]] <- cbind(var.info, n.obs, freq, weight)
                  
                  if (n.site > 0) {
                      # mean impute missing values
                      if (any(n.obs < nrow(geno))) {
                          geno <- .meanImpute(geno, freq)
                      }

                      # do the test
                      assoc <- testVariantSet(null.model, G=geno, weights=weight, test=test, 
                                              burden.test=burden.test, rho=rho,
                                              pval.method=pval.method)
                      res[[i]] <- cbind(res[[i]], assoc)
                  }

                  if (verbose & i %% 100 == 0) {
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
                   test=c("Burden", "SKAT", "SMMAT"),
                   burden.test=c("Score", "Wald"), rho=0,
                   pval.method=c("davies", "kuonen", "liu"),
                   verbose=TRUE) {

              # check argument values
              test <- match.arg(test)
              burden.test <- match.arg(burden.test)
              pval.method <- match.arg(pval.method)

              # filter samples to match null model
              sample.id <- null.model$sample.id
              if (!is.null(sample.id)) {
                  sample.index <- match(sample.id, getScanID(gdsobj))
              } else {
                  sample.index <- match(rownames(null.model$model.matrix),
                                        sampleNames(getScanAnnotation(gdsobj)))
                  sample.id <- getScanID(gdsobj)[sample.index]
              }

              # results
              res <- list()
              res.var <- list()
              i <- 1
              n.iter <- length(snpFilter(gdsobj))
              iterate <- TRUE
              while (iterate) {
                  var.info <- data.frame(variant.id=getSnpID(gdsobj),
                                         chr=getChromosome(gdsobj, char=TRUE),
                                         pos=getPosition(gdsobj),
                                         stringsAsFactors=FALSE)
                  
                  geno <- getGenotypeSelection(gdsobj, scanID=sample.id, order="selection",
                                               transpose=TRUE)
                  
                  # allele frequency
                  freq <- .alleleFreq(gdsobj, geno, sample.index=sample.index)
                  # exclude monomorphic variants
                  mono <- freq %in% c(0,1)
                  # exclude variants with freq > max
                  excl <-  mono | freq > AF.max | is.na(freq)
                  if (any(excl)) {
                      var.info <- var.info[!excl,,drop=FALSE]
                      geno <- geno[,!excl,drop=FALSE]
                      freq <- freq[!excl]
                  }

                  # number of variant sites
                  n.site <- length(unique(var.info$variant.id))

                  # number of alternate alleles
                  n.alt <- sum(geno, na.rm=TRUE)
                  
                  # number of samples with observed alternate alleles > 0
                  n.sample.alt <- sum(rowSums(geno, na.rm=TRUE) > 0)
               
                  # number of non-missing samples
                  n.obs <- colSums(!is.na(geno))
                  
                  # weights
                  if (is.null(weight.user)) {
                      # Beta weights
                      weight <- .weightFromFreq(freq, weight.beta)
                  } else {
                      # user supplied weights
                      weight <- getSnpVariable(gdsobj, weight.user)
                      weight <- weight[!excl]
                  }
                  
                  res[[i]] <- data.frame(n.site, n.alt, n.sample.alt)
                  res.var[[i]] <- cbind(var.info, n.obs, freq, weight)
                  
                  if (n.site > 0) {
                      # mean impute missing values
                      if (any(n.obs < nrow(geno))) {
                          geno <- .meanImpute(geno, freq)
                      }

                      # do the test
                      assoc <- testVariantSet(null.model, G=geno, weights=weight, test=test, 
                                              burden.test=burden.test, rho=rho,
                                              pval.method=pval.method)
                      res[[i]] <- cbind(res[[i]], assoc)
                  }

                  if (verbose & i %% 100 == 0) {
                      message(paste("Iteration", i , "of", n.iter, "completed"))
                  }
                  i <- i + 1
                  iterate <- GWASTools::iterateFilter(gdsobj)
              }

              res <- list(results=bind_rows(res), variantInfo=res.var)
              .annotateAssoc(gdsobj, res)
          })
