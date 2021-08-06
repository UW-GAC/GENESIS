setGeneric("assocTestSingle", function(gdsobj, ...) standardGeneric("assocTestSingle"))

## do we want the GxE.return.cov option?
## do we want to make imputing to the mean optional?
setMethod("assocTestSingle",
          "SeqVarIterator",
          function(gdsobj, null.model, test=c("Score", "Score.SPA", "BinomiRare", "CMP"),
                   recalc.pval.thresh=0.05, fast.score.SE=FALSE, GxE=NULL,
                   sparse=TRUE, imputed=FALSE, male.diploid=TRUE, genome.build=c("hg19", "hg38"), 
                   BPPARAM=bpparam(), verbose=TRUE) {
              test <- match.arg(test)

              # don't use sparse matrices for imputed dosages
              if (imputed) sparse <- FALSE

              # Convert old null model format if necessary.
              null.model <- .updateNullModelFormat(null.model)

              # check that the provided null model is compatible with the requested test
              .checkNullModelTestSingle(null.model = null.model, test = test, 
              	recalc.pval.thresh = recalc.pval.thresh, fast.score.SE = fast.score.SE, GxE = GxE)
              
              # coerce null.model if necessary
              if (sparse) null.model <- .nullModelAsMatrix(null.model)

              # filter samples to match null model
              sample.index <- .setFilterNullModel(gdsobj, null.model, verbose=verbose)
              
              # get sex for calculating allele freq
              sex <- validateSex(gdsobj)[sample.index]
              
              if (!is.null(GxE)) GxE <- .modelMatrixColumns(null.model, GxE)
              
              # check ploidy
              if (SeqVarTools:::.ploidy(gdsobj) == 1) male.diploid <- FALSE

              # results
              #n.iter <- length(variantFilter(gdsobj))
              #set.messages <- ceiling(n.iter / 100) # max messages = 100
                  
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
                  
                  var.info <- variantInfo(gdsobj, alleles=FALSE, expanded=TRUE)
                  
                  if (!imputed) {
                      geno <- expandedAltDosage(gdsobj, use.names=FALSE, sparse=sparse)[sample.index,,drop=FALSE]
                  } else {
                      geno <- imputedDosage(gdsobj, use.names=FALSE)[sample.index,,drop=FALSE]
                  }
                  
                  chr <- chromWithPAR(gdsobj, genome.build=genome.build)
                  
                  return(list(var.info=var.info, geno=geno, chr=chr))
              }

              res <- bpiterate(ITER, .testGenoBlockSingle, BPPARAM=BPPARAM,
                               sex=sex, null.model=null.model, test=test,
                               recalc.pval.thresh=recalc.pval.thresh, 
                               fast.score.SE=fast.score.SE, GxE=GxE,
                               sparse=sparse, imputed=imputed,
                               male.diploid=male.diploid)
              .stopOnError(res)
              as.data.frame(rbindlist(res))
          })



setMethod("assocTestSingle",
          "GenotypeIterator",
          function(gdsobj, null.model, test=c("Score", "Score.SPA", "BinomiRare", "CMP"),
                   recalc.pval.thresh=0.05, GxE=NULL, 
                   male.diploid=TRUE, BPPARAM=bpparam(), verbose=TRUE) {
              test <- match.arg(test)

              # Convert old null model format if necessary.
              null.model <- .updateNullModelFormat(null.model)

              # check that the provided null model is compatible with the requested test
              .checkNullModelTestSingle(null.model = null.model, test = test, 
              	recalc.pval.thresh = recalc.pval.thresh, fast.score.SE = FALSE, GxE = GxE)

              # filter samples to match null model
              sample.index <- .sampleIndexNullModel(gdsobj, null.model)
              
              # get sex for calculating allele freq
              sex <- validateSex(gdsobj)[sample.index]
              
              if (!is.null(GxE)) GxE <- .modelMatrixColumns(null.model, GxE)

              # results
              # n.iter <- length(snpFilter(gdsobj))
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
                  
                  return(list(var.info=var.info, geno=geno, chr=var.info$chr))
              }
              
              res <- bpiterate(ITER, .testGenoBlockSingle, BPPARAM=BPPARAM,
                               sex=sex, null.model=null.model, test=test,
                               recalc.pval.thresh=recalc.pval.thresh, 
                               fast.score.SE=FALSE, GxE=GxE,
                               sparse=FALSE, imputed=FALSE,
                               male.diploid=male.diploid)
              .stopOnError(res)
              as.data.frame(rbindlist(res))
          })



# function to process a block of genotype data
.testGenoBlockSingle <- function(x, sex, null.model, test,
                                 recalc.pval.thresh, fast.score.SE, GxE,
                                 sparse, imputed, male.diploid, ...) {

    x <- .prepGenoBlock(x, sex=sex, test=test, imputed=imputed, male.diploid=male.diploid)
    var.info <- x$var.info
    n.obs <- x$n.obs
    freq <- x$freq
    geno <- x$geno
    rm(x)
    
    # do the test
    if (ncol(geno) == 0){
        res.i <- NULL
    } else {
        assoc <- testGenoSingleVar(null.model, G=geno, E=GxE, test=test,
                                   recalc.pval.thresh=recalc.pval.thresh,
                                   fast.score.SE=fast.score.SE)
        
        res.i <- cbind(var.info, n.obs, freq, assoc)
    }
    
    #if (verbose & n.iter > 1 & i %% set.messages == 0) {
    #    message(paste("Iteration", i , "of", n.iter, "completed"))
    #}
    return(res.i)
}


# function to pre-process genotype data before testing
.prepGenoBlock <- function(x, sex, test, imputed, male.diploid) {
    
    var.info <- x$var.info
    geno <- x$geno
    chr <- x$chr
    rm(x)
    
    # take note of number of non-missing samples
    #n.obs <- colSums(!is.na(geno))
    n.obs <- .countNonMissing(geno, MARGIN = 2)
    
    # allele frequency
    freq <- .alleleFreq(geno, chr, sex, male.diploid=male.diploid)
    
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
    
    # mean impute missing values
    if (any(n.obs < nrow(geno))) {
        geno <- .meanImpute(geno, freq$freq)
    }
    
    return(list(var.info=var.info, n.obs=n.obs, freq=freq, geno=geno))
}

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
