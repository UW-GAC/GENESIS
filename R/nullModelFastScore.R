## fit the null model
## estimate the ratio (se.ratio) of the true SE(score) to the fast SE(score)
## update the null model to contain this parameter
## this needs to be run before using the fast.score.se option in assocTestSingle
## this is based off of the SAIGE variance approximation
setGeneric("fitNullModelFastScore", function(x, ...) standardGeneric("fitNullModelFastScore"))

setMethod("fitNullModelFastScore",
          "SeqVarData",
          function(x, outcome,
                   covars = NULL,
                   cov.mat = NULL,
                   group.var = NULL,
                   family = "gaussian",
                   two.stage = FALSE,
                   norm.option = c("all", "by.group"),
                   rescale = c("none", "model", "residSD"),
                   start = NULL,
                   AIREML.tol = 1e-4,
                   max.iter = 100,
                   EM.iter = 0,
                   drop.zeros = TRUE,
                   return.small = TRUE,
                   variant.id = NULL,
                   nvar = 100, 
                   min.mac = 20, 
                   sparse = TRUE, 
                   imputed = FALSE, 
                   male.diploid = TRUE, 
                   genome.build = c("hg19", "hg38"), 
                   verbose=TRUE) {

              # check if this is a mixed model
              if(is.null(cov.mat)){
                  stop('Fast Score.SE approximation only applies to mixed models (cov.mat not NULL); use fitNullModel for the specified model')
              }

              # fit the null model
              null.model <- fitNullModel(sampleData(x), outcome = outcome, covars = covars, cov.mat = cov.mat,
                                         group.var = group.var, family = family, two.stage = two.stage,
                                         norm.option = norm.option, rescale = rescale, start = start,
                                         AIREML.tol = AIREML.tol, max.iter = max.iter, EM.iter = EM.iter,
                                         drop.zeros = drop.zeros, return.small = FALSE, verbose = verbose)

              # calculate true score SE and the fast approximation
              tab <- calcScore(gdsobj = x, null.model = null.model,
                               variant.id = variant.id, nvar = nvar, min.mac = min.mac,
                               sparse = sparse, imputed = imputed, male.diploid = male.diploid,
                               genome.build = genome.build, verbose = verbose)

              # update the null model with the se.correction factor
              null.model <- nullModelFastScore(null.model = null.model, tab = tab, return.small = return.small)

              null.model
          })


## calculate the Score and Score.SE for a specified or random set of variants
## when a mixed model, also calculates the fast Score.SE approximation and the ratio to the true Score.SE
calcScore <- function(gdsobj,
                      null.model, 
                      variant.id = NULL,
                      nvar = 100, 
                      min.mac = 20, 
                      sparse=TRUE, 
                      imputed=FALSE, 
                      male.diploid=TRUE, 
                      genome.build=c("hg19", "hg38"), 
                      verbose=TRUE){

     # Update null model format
     null.model <- .updateNullModelFormat(null.model)

     # samples in null model
     sampid <- null.model$fit$sample.id
     if(verbose) message(paste('null.model has', length(sampid), 'samples'))

     if(is.null(variant.id)){
          # select a random set of variants meeting min.mac in the sample set
          variant.id <- .selectRandomVars(gdsobj, sample.id = sampid, nvar = nvar, min.mac = min.mac, verbose = verbose)
          if(verbose) message(paste('Selected', length(variant.id), 'random variants with MAC >=', min.mac))
     }else{
          if(verbose) message(paste('User provided', length(variant.id), 'variants in variant.id'))
     }
     
     # calculate Score.SE using both approaches
     tab <- .calcScore(gdsobj, null.model, variant.id = variant.id, sparse = sparse, imputed = imputed, 
                         male.diploid = male.diploid, genome.build = genome.build, verbose = verbose)
     return(tab)
}


## updates the null model object with the parameter needed for fast.score.se	
nullModelFastScore <- function(null.model, tab, return.small = TRUE){
    # Update null model format
     null.model <- .updateNullModelFormat(null.model)

    # rbind a list of tables
    if(class(tab) == 'list') tab <- data.table::rbindlist(tab)

    # compute the mean of se.ratio
    r <- mean(tab$se.ratio)

    # compute SE(r) and check if enough variants were included
    r.se <- sd(tab$se.ratio)/sqrt(nrow(tab))
    val <- r.se/r
    message(paste('mean se.ratio: r = ', r))
    message(paste('SE(r)/r = ', val))
    if(val > 0.0025){
        warning(paste('It is recommended that SE(r)/r be < 0.0025. It is suggested to increase the number of variants in tab; try at least', 
          round(nrow(tab)*(val/0.0025)^2, 0), 'variants'))
    }

    # update the null model
    if(return.small){
        null.model <- nullModelSmall(null.model)
    }
    null.model$se.correction <- r
    null.model$se.table <- tab
    return(null.model)
}


## function to get a random sample of variants with a minimum MAC
## this is called by calcScore
setGeneric(".selectRandomVars", function(gdsobj, ...) standardGeneric(".selectRandomVars"))

setMethod(".selectRandomVars",
          "SeqVarData",
          function(gdsobj, 
                   sample.id = NULL, 
                   nvar = 100, 
                   min.mac = 20, 
                   verbose = TRUE){

            # filter to sample set in null.model
            SeqArray::seqSetFilter(gdsobj, sample.id = sample.id, verbose = verbose)

            # get the initial set of filtered variants 
            # (accounts for any pre-filtering; e.g. by PASS status)
            var.filt <- which(SeqArray::seqGetFilter(gdsobj)$variant.sel)

            # number of variants in the GDS object
            nvar.gds <- length(var.filt)
            if(nvar.gds < nvar) stop('requested more variants than available in gdsobj')

            out <- NULL
            while(length(out) < nvar){
              # sample variants
              var.rand <- sort(sample(var.filt, size = nvar, replace = FALSE))
              # filter to random sample
              SeqArray::seqSetFilter(gdsobj, variant.sel = var.rand, verbose = verbose)
              # filter to MAC threshold
              if(min.mac > 0) SeqArray::seqSetFilterCond(gdsobj, mac = min.mac, verbose = verbose)
              # collect selected variants
              out <- sort(unique(c(out, seqGetData(gdsobj, 'variant.id'))))
            }

            # sample down to number requested
            out <- sample(out, size = nvar, replace = FALSE)
            return(out)
          })


## function to calculate both the true and fast Score.SE for a specified set of variants
## computes the ratio of the two Score.SE estimates for calculating the correction factor
## this is called by calcScore
setGeneric(".calcScore", function(gdsobj, ...) standardGeneric(".calcScore"))

setMethod(".calcScore",
          "SeqVarData",
          function(gdsobj, 
                   null.model, 
                   variant.id,
                   sparse=TRUE, 
                   imputed=FALSE, 
                   male.diploid=TRUE, 
                   genome.build=c("hg19", "hg38"), 
                   verbose = TRUE){

               # don't use sparse matrices for imputed dosages
               if (imputed) sparse <- FALSE

               # coerce null.model if necessary
               if (sparse) null.model <- .nullModelAsMatrix(null.model)

               # filter samples to match null model
               sample.index <- .setFilterNullModel(gdsobj, null.model, verbose=verbose)

               # filter variants
               seqSetFilter(gdsobj, variant.id = variant.id, verbose = verbose)

               # check ploidy
               if (SeqVarTools:::.ploidy(gdsobj) == 1) male.diploid <- FALSE

               # variant info
               var.info <- variantInfo(gdsobj, alleles=FALSE, expanded=TRUE)

               # read in genotype data
               if (!imputed) {
                    geno <- expandedAltDosage(gdsobj, use.names=FALSE, sparse=sparse)[sample.index,,drop=FALSE]
               } else {
                    geno <- imputedDosage(gdsobj, use.names=FALSE)[sample.index,,drop=FALSE]
               }

               # take note of number of non-missing samples
               n.obs <- .countNonMissing(geno, MARGIN = 2)

               # allele frequency
               freq <- .alleleFreq(gdsobj, geno, sample.index=sample.index,
                                        male.diploid=male.diploid, genome.build=genome.build)

               # filter monomorphic variants
               keep <- .filterMonomorphic(geno, count=n.obs, freq=freq$freq, imputed=imputed)
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

               geno <- .genoAsMatrix(null.model, geno)

               # score
               score <- as.vector(crossprod(geno, null.model$fit$resid.PY))

               # true variance
               Gtilde <- calcGtilde(null.model, geno)
               se.true <- sqrt(colSums(Gtilde^2)) # sqrt(GPG)

               # results
               tab <- cbind(var.info, n.obs, freq, Score = score, Score.SE = se.true)

               if(null.model$model$family$mixedmodel){
                    # approx SE
                    Gtilde <- calcGtildeFast(null.model, geno, r = 1)
                    se.fast <- sqrt(colSums(Gtilde^2)) #sqrt(GWG)

                    # results
                    tab <- cbind(tab, Score.SE.fast = se.fast, se.ratio = se.true/se.fast)
               }
               
               return(tab)
          })

