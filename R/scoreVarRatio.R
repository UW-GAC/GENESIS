## function that estimates the ratio (r) of the true Var(score) to the approximated Var(score)
## this needs to be run before using the approx.score.var option in assocTestSingle
## this is based off of the SAIGE variance approximation

setGeneric("scoreVarRatio", function(gdsobj, ...) standardGeneric("scoreVarRatio"))

setMethod("scoreVarRatio",
          "SeqVarIterator",
          function(gdsobj, null.model, sparse=TRUE, imputed=FALSE, 
          			male.diploid=TRUE, genome.build=c("hg19", "hg38"), verbose=TRUE){

          	# don't use sparse matrices for imputed dosages
          	if (imputed) sparse <- FALSE

               # Convert old null model format if necessary.
               null.model <- .updateNullModelFormat(null.model)

          	# coerce null.model if necessary
          	if (sparse) null.model <- .nullModelAsMatrix(null.model)

          	# filter samples to match null model
          	sample.index <- .setFilterNullModel(gdsobj, null.model, verbose=verbose)

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

          	geno <- .genoAsMatrix(nullmod, geno)

          	# true variance
          	Gtilde <- calcGtilde(nullmod, geno)
          	GPG <- colSums(Gtilde^2)

          	# approx variance
          	Gtilde <- calcGtildeApprox(nullmod, geno, r = 1)
          	GWG <- colSums(Gtilde^2)

          	# results
          	res <- cbind(var.info, n.obs, freq, var.true = GPG, var.approx = GWG, r = GPG/GWG)
          	r <- mean(res$r)

          })

setMethod("scoreVarRatio",
          "GenotypeIterator",
          function(gdsobj, null.model){

          })