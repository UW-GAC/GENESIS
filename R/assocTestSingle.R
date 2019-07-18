setGeneric("assocTestSingle", function(gdsobj, ...) standardGeneric("assocTestSingle"))

## do we want the GxE.return.cov option?
## do we want to make imputing to the mean optional?
setMethod("assocTestSingle",
          "SeqVarIterator",
          function(gdsobj, null.model, test=c("Score", "Wald"), GxE=NULL, sparse=TRUE, imputed=FALSE, verbose=TRUE) {
              test <- match.arg(test)

              # don't use sparse matrices for imputed dosages
              if (imputed) sparse <- FALSE

              # coerce null.model if necessary
              if (sparse) null.model <- .nullModelAsMatrix(null.model)
              
              # filter samples to match null model
              sample.index <- .setFilterNullModel(gdsobj, null.model, verbose=verbose)

              if (!is.null(GxE)) GxE <- .modelMatrixColumns(null.model, GxE)
              
              # results
              res <- list()
              n.iter <- length(variantFilter(gdsobj))
              set.messages <- ceiling(n.iter / 100) # max messages = 100
              i <- 1
              iterate <- TRUE
              while (iterate) {
                  var.info <- variantInfo(gdsobj, alleles=FALSE, expanded=TRUE)

                  if (!imputed) {
                      geno <- expandedAltDosage(gdsobj, use.names=FALSE, sparse=sparse)[sample.index,,drop=FALSE]
                  } else {
                      geno <- imputedDosage(gdsobj, use.names=FALSE)[sample.index,,drop=FALSE]
                  }
                  
                  # take note of number of non-missing samples
                  n.obs <- colSums(!is.na(geno))
                  
                  # allele frequency
                  freq <- .alleleFreq(gdsobj, geno, sample.index=sample.index)
                  
                  # filter monomorphic variants
                  keep <- .filterMonomorphic(geno, count=n.obs, freq=freq, imputed=imputed)
                  if (!all(keep)) {
                      var.info <- var.info[keep,,drop=FALSE]
                      geno <- geno[,keep,drop=FALSE]
                      n.obs <- n.obs[keep]
                      freq <- freq[keep]
                  }

                  # mean impute missing values
                  if (any(n.obs < nrow(geno))) {
                      geno <- .meanImpute(geno, freq)
                  }

                  # do the test
                  assoc <- testGenoSingleVar(null.model, G=geno, E=GxE, test=test)

                  res[[i]] <- cbind(var.info, n.obs, freq, assoc)
                  
                  if (verbose & n.iter > 1 & i %% set.messages == 0) {
                      message(paste("Iteration", i , "of", n.iter, "completed"))
                  }
                  i <- i + 1
                  iterate <- iterateFilter(gdsobj, verbose=FALSE)
              }

              do.call(rbind, res)
          })



setMethod("assocTestSingle",
          "GenotypeIterator",
          function(gdsobj, null.model, test=c("Score", "Wald"), GxE=NULL, verbose=TRUE) {
              test <- match.arg(test)

              # filter samples to match null model
              sample.index <- .sampleIndexNullModel(gdsobj, null.model)
              
              if (!is.null(GxE)) GxE <- .modelMatrixColumns(null.model, GxE)
              
              # results
              res <- list()
              n.iter <- length(snpFilter(gdsobj))
              set.messages <- ceiling(n.iter / 100) # max messages = 100
              i <- 1
              iterate <- TRUE
              while (iterate) {
                  var.info <- variantInfo(gdsobj)
                  
                  geno <- getGenotypeSelection(gdsobj, scan=sample.index, order="selection",
                                               transpose=TRUE, use.names=FALSE, drop=FALSE)
                  
                  # take note of number of non-missing samples
                  n.obs <- colSums(!is.na(geno))
                  
                  # allele frequency
                  freq <- .alleleFreq(gdsobj, geno, sample.index=sample.index)
                  
                  # filter monomorphic variants
                  keep <- .filterMonomorphic(geno, count=n.obs, freq=freq)
                  if (!all(keep)) {
                      var.info <- var.info[keep,,drop=FALSE]
                      geno <- geno[,keep,drop=FALSE]
                      n.obs <- n.obs[keep]
                      freq <- freq[keep]
                  }

                  # mean impute missing values
                  if (any(n.obs < nrow(geno))) {
                      geno <- .meanImpute(geno, freq)
                  }

                  # do the test
                  assoc <- testGenoSingleVar(null.model, G=geno, E=GxE, test=test)

                  res[[i]] <- cbind(var.info, n.obs, freq, assoc)
                  
                  if (verbose & n.iter > 1 & i %% set.messages == 0) {
                      message(paste("Iteration", i , "of", n.iter, "completed"))
                  }
                  i <- i + 1
                  iterate <- GWASTools::iterateFilter(gdsobj)
              }

              do.call(rbind, res)
          })
