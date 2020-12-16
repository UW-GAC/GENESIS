## estimate the ratio (score.se.ratio) of the true SE(score) to the fast SE(score)
## update the null model to contain this parameter
## this needs to be run before using the approx.score.se option in assocTestSingle
## this is based off of the SAIGE variance approximation

setGeneric("updateNullModelApproxScoreSE", function(gdsobj, ...) standardGeneric("updateNullModelApproxScoreSE"))

setMethod("updateNullModelApproxScoreSE",
          "SeqVarData",
          function(gdsobj, null.model, nvar = 100, min.mac = 20, 
                    sparse=TRUE, imputed=FALSE, male.diploid=TRUE, genome.build=c("hg19", "hg38"), 
                    verbose=TRUE){

               # check for W matrix
               if (is.null(null.model$W)) stop('This null model was created with an older version of GENESIS and is not compatible with this analysis; 
                    please re-run your null model with the latest version.')

               # samples in null model
               sampid <- null.model$fit$sample.id
               if(verbose) message(paste('null.model has', length(sampid), 'samples'))

               # select a random set of variants meeting min.mac in the sample set
               varid <- .selectRandomVars(gdsobj, sample.id = sampid, nvar = nvar, min.mac = min.mac)
               if(verbose) message(paste('Selected', length(varid), 'variants with MAC >=', min.mac))

               # calculate Score.SE using both approaches
               tab <- .scoreSEtable(gdsobj, null.model, variant.id = varid, 
                                   sparse = sparse, imputed = imputed, male.diploid = male.diploid, genome.build = genome.build, verbose = verbose)
          	
               # calculate the score.se.ratio estimate
               score.se.ratio <- .scoreSEratio(tab)

               # update the null model
               null.model <- nullModelSmall(null.model)
               null.model$score.se.ratio <- score.se.ratio
               return(null.model)
          })

setMethod("updateNullModelApproxScoreSE",
          "GenotypeIterator",
          function(gdsobj){

          })



## function to get a random sample of variants with a minimum MAC
## this is called by updateNullModelApproxScoreSE
setGeneric(".selectRandomVars", function(gdsobj, ...) standardGeneric(".selectRandomVars"))

setMethod(".selectRandomVars",
          "SeqVarData",
          function(gdsobj, sample.id = NULL, nvar = 100, min.mac = 20){

               # number of variants in the GDS object
               nvar.gdsobj <- SeqArray::seqSummary(gdsobj, verbose = FALSE)$num.variant
               if(nvar.gdsobj < nvar) stop('requested more variants than available in gdsobj')

               out <- NULL
               while(length(out) < nvar){
                    # reset filters
                    SeqArray::seqSetFilter(gdsobj, sample.id = sample.id, verbose = FALSE)
                    # sample variants
                    var.rand <- sort(sample.int(nvar.gdsobj, size = nvar, replace = FALSE))
                    # filter to random sample
                    SeqArray::seqSetFilter(gdsobj, variant.sel = var.rand, verbose = FALSE)
                    # filter to MAC threshold
                    if(min.mac > 0) SeqArray::seqSetFilterCond(gdsobj, mac = min.mac, verbose = FALSE)
                    # collect selected variants
                    out <- sort(unique(c(out, seqGetData(gdsobj, 'variant.id'))))
               }

               # sample down to number requested
               out <- sample(out, size = nvar, replace = FALSE)
               return(out)
               })

setMethod(".selectRandomVars",
          "GenotypeIterator",
          function(gdsobj){

          })


## function to calculate both the true and fast Score.SE for a specified set of variants
## computes the ratio of the two Score.SE estimates for calculating the correction factor
## this is called by updateNullModelApproxScoreSE
setGeneric(".scoreSEtable", function(gdsobj, ...) standardGeneric(".scoreSEtable"))

setMethod(".scoreSEtable",
          "SeqVarData",
          function(gdsobj, null.model, variant.id,
                    sparse=TRUE, imputed=FALSE, male.diploid=TRUE, genome.build=c("hg19", "hg38"), verbose = TRUE){

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

               # true variance
               Gtilde <- calcGtilde(null.model, geno)
               se.true <- sqrt(colSums(Gtilde^2)) # sqrt(GPG)

               # approx variance
               Gtilde <- calcGtildeApprox(null.model, geno, r = 1)
               se.approx <- sqrt(colSums(Gtilde^2)) #sqrt(GWG)

               # results
               tab <- cbind(var.info, n.obs, freq, score.SE.true = se.true, score.SE.approx = se.approx, score.SE.ratio = se.true/se.approx)
               return(tab)

          })

## computes the score.se.ratio estimate across all variants in the variant table
.scoreSEratio <- function(x){
     # compute the mean of score.se.ratio
     mu <- mean(x$score.SE.ratio)

     # compute SE(mu) and check if enough variants were included
     mu.se <- sd(x$score.SE.ratio)/sqrt(nrow(x))
     val <- mu.se/mu
     message(paste('mean score.SE.ratio: mu = ', mu))
     message(paste('SE(mu)/mu = ', val))
     if(val > 0.0025) warning(paste('It is recommended that SE(mu)/mu be < 0.0025. It is suggested to increase the number of variants in tab; try at least', 
          round(nrow(x)*(val/0.0025)^2, 0), 'variants'))

     return(mu)
}


