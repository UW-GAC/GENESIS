setGeneric("assocTestSingle", function(gdsobj, ...) standardGeneric("assocTestSingle"))

## do we want the GxE.return.cov option?
## do we want to make imputing to the mean optional?
setMethod("assocTestSingle",
          "SeqVarIterator",
          function(gdsobj, null.model, test=c("Score", "Score.SPA", "BinomiRare", "CMP"),
                   recalc.pval.thresh=0.05, fast.score.SE=FALSE, GxE=NULL,
                   sparse=TRUE, imputed=FALSE, male.diploid=TRUE, genome.build=c("hg19", "hg38"), 
                   verbose=TRUE, geno.test=c('additive','dominant','recessive')) {

              test <- match.arg(test)

              # don't use sparse matrices for imputed dosages
              if (imputed) sparse <- FALSE

              # Convert old null model format if necessary.
              null.model <- .updateNullModelFormat(null.model)

              # check that the provided null model is compatible with the requested test
              .checkNullModelTestSingle(null.model = null.model, test = test, 
                                        recalc.pval.thresh = recalc.pval.thresh, fast.score.SE = fast.score.SE, GxE = GxE)
              ## could add a recessive + impute, etc. checking to this function
              
              # coerce null.model if necessary
              if (sparse) null.model <- .nullModelAsMatrix(null.model)

              # filter samples to match null model
              sample.index <- .setFilterNullModel(gdsobj, null.model, verbose=verbose)
              if (!is.null(GxE)) GxE <- .modelMatrixColumns(null.model, GxE)

              # check ploidy
              if (SeqVarTools:::.ploidy(gdsobj) == 1) male.diploid <- FALSE

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
                  #n.obs <- colSums(!is.na(geno))
                  n.obs <- .countNonMissing(geno, MARGIN = 2)

                  if(geno.test=='recessive'){ # is recessive
                      ## if wanting to test a recessive model, the genotypes 0 and 1 are '0' and the genotype 2 is '1',
                      ##eg. indicator of participant having 2 copies of the minor allele
                      geno[geno==1] <- 0
                      geno[geno==2] <- 1
                  }else if(geno.test=='dominant'){
                      geno[geno==2] <- 1
                  }
                  
                  # allele frequency
                  freq <- .alleleFreq(gdsobj, geno, sample.index=sample.index,
                                      male.diploid=male.diploid, genome.build=genome.build, geno.test=geno.test)
                  
                  # filter monomorphic variants
                  keep <- .filterMonomorphic(geno, count=n.obs, freq=freq$freq, imputed=imputed)

                  # for BinomiRare and CMP, restrict to variants where the alternate allele is minor
                  if (test %in% c("BinomiRare", "CMP")) {
                      keep <- keep & (freq$freq <= 0.5)
                  }

                  if (!all(keep)) {
                      var.info <- var.info[keep,,drop=FALSE]
                      geno <- geno[,keep,drop=FALSE]
                      n.obs <- n.obs[keep]
                      freq <- freq[keep,,drop=FALSE]
                  }

                  ## mean impute missing values
                  if (any(n.obs < nrow(geno))) {
                      geno <- .meanImpute(geno, freq$freq)
                  }

                  ## do the test
                  if (ncol(geno) == 0){
                      res[[i]] <- NULL
                  } else {
                      assoc <- testGenoSingleVar(null.model, G=geno, E=GxE, test=test,
                                                 recalc.pval.thresh=recalc.pval.thresh,
                                                 fast.score.SE=fast.score.SE)
                      ## description of the freq and MAC columns will need to change under the recessive coding options

                      res[[i]] <- cbind(var.info, n.obs, freq, assoc)
                  }

                  if (verbose & n.iter > 1 & i %% set.messages == 0) {
                      message(paste("Iteration", i , "of", n.iter, "completed"))
                  }
                  message('end of first iteration loop')
                  i <- i + 1
                  iterate <- iterateFilter(gdsobj, verbose=FALSE)
              }

              as.data.frame(rbindlist(res))
          })



setMethod("assocTestSingle",
          "GenotypeIterator",
          function(gdsobj, null.model, test=c("Score", "Score.SPA", "BinomiRare", "CMP"),
                   recalc.pval.thresh=0.05, GxE=NULL, 
                   male.diploid=TRUE, verbose=TRUE) {
              test <- match.arg(test)

              # Convert old null model format if necessary.
              null.model <- .updateNullModelFormat(null.model)

              # check that the provided null model is compatible with the requested test
              .checkNullModelTestSingle(null.model = null.model, test = test, 
              	recalc.pval.thresh = recalc.pval.thresh, fast.score.SE = FALSE, GxE = GxE)

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
                  #n.obs <- colSums(!is.na(geno))
                  n.obs <- .countNonMissing(geno, MARGIN = 2)

                  # allele frequency
                  freq <- .alleleFreq(gdsobj, geno, sample.index=sample.index,
                                      male.diploid=male.diploid)

                  # filter monomorphic variants (and max alternate frequency variants)
                  keep <- .filterMonomorphic(geno, count=n.obs, freq=freq$freq)

                  # for BinomiRare and CMP, restrict to variants where the alternate allele is minor
                  if (test %in% c("BinomiRare", "CMP")) {
                      keep <- keep & (freq$freq <= 0.5)
                  }

                  if (!all(keep)) {
                      var.info <- var.info[keep,,drop=FALSE]
                      geno <- geno[,keep,drop=FALSE]
                      n.obs <- n.obs[keep]
                      freq <- freq[keep,,drop=FALSE]
                  }

                  # mean impute missing values
                  if (any(n.obs < nrow(geno))) {
                      geno <- .meanImpute(geno, freq$freq)
                  }

                  # do the test
                  if (ncol(geno) == 0){
                      res[[i]] <- NULL
                  } else {
                      assoc <- testGenoSingleVar(null.model, G=geno, E=GxE, test=test,
                                                 recalc.pval.thresh=recalc.pval.thresh)

                      res[[i]] <- cbind(var.info, n.obs, freq, assoc)
                  }

                  if (verbose & n.iter > 1 & i %% set.messages == 0) {
                      message(paste("Iteration", i , "of", n.iter, "completed"))
                  }
                  i <- i + 1
                  iterate <- GWASTools::iterateFilter(gdsobj)
              }

              as.data.frame(rbindlist(res))
          })


# check that the provided null model is compatible with the requested test
.checkNullModelTestSingle <- function(null.model, test, recalc.pval.thresh, fast.score.SE, GxE){
	calc.score <- test %in% c("Score", "Score.SPA") | (recalc.pval.thresh < 1)

	if(fast.score.SE && !isNullModelFastScore(null.model)){
		stop("null.model must have se.correction when fast.score.SE = TRUE; re-fit your null.model using `fitNullModelFastScore` or update your null.model using `nullModelFastScore`")
	}

	if(calc.score && !(fast.score.SE) && isNullModelSmall(null.model)){
		stop("small null.model cannot be used with test options provided")
	}

	if(!is.null(GxE) && isNullModelSmall(null.model)){
		stop("small null.model cannot be used with GxE")
	}
}
