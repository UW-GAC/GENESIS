setGeneric("assocTestSingle", function(gdsobj, ...) standardGeneric("assocTestSingle"))

## do we want the GxE.return.cov option?
## do we want to make imputing to the mean optional?
setMethod("assocTestSingle",
          "SeqVarIterator",
          function(gdsobj, null.model, test=c("Score", "Wald"), GxE=NULL, verbose=TRUE) {
              test <- match.arg(test)

              # filter samples to match null model
              sample.index <- .setFilterNullModel(gdsobj, null.model, verbose=verbose)
              
              # results
              res <- list()
              n.iter <- length(variantFilter(gdsobj))
              i <- 1
              iterate <- TRUE
              while (iterate) {
                  var.info <- variantInfo(gdsobj, alleles=FALSE, expanded=TRUE)
                  
                  geno <- expandedAltDosage(gdsobj, use.names=FALSE, sparse=TRUE)[sample.index,,drop=FALSE]
                  
                  # allele frequency
                  freq <- .alleleFreq(gdsobj, geno)

                  # take note of number of non-missing samples
                  n.obs <- colSums(!is.na(geno))
                  
                  # mean impute missing values
                  if (any(n.obs < nrow(geno))) {
                      geno <- .meanImpute(geno, freq)
                  }

                  # do the test
                  if (!is.null(GxE)) GxE <- .modelMatrixColumns(null.model, GxE)
                  assoc <- testGenoSingleVar(null.model, G=geno, E=GxE, test=test)
                  # set monomorphs to NA - do we want to skip testing these to save time?
                  assoc[freq %in% c(0,1),] <- NA

                  res[[i]] <- cbind(var.info, n.obs, freq, assoc)
                  
                  if (verbose & i %% 100 == 0) {
                      message(paste("Iteration", i , "of", n.iter, "completed"))
                  }
                  i <- i + 1
                  iterate <- iterateFilter(gdsobj, verbose=FALSE)
              }

              do.call(rbind, res)
          })
