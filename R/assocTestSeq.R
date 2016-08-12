assocTestSeq <- function(seqData,
                         nullModObj,
                         aggVarList,
                         AF.sample = NULL,
                         AF.range = c(0,1),
                         weight.beta = c(0.5,0.5),
                         weight.user = NULL,
                         test = "Burden",
                         burden.test = "Score",
                         rho = 0,
                         pval.method = "kuonen",
                         verbose = TRUE){

	# save the filter
	seqFilt.original <- seqGetFilter(seqData)
        # reset so indexing works
        seqResetFilter(seqData, verbose=FALSE)

	# check the parameters
	param <- .paramChecks(seqData = seqData, AF.range = AF.range, weight.beta = weight.beta, weight.user = weight.user, 
                              test = test, burden.test = burden.test, family = nullModObj$family$family, 
                              mixedmodel = nullModObj$family$mixedmodel, rho = rho, pval.method = pval.method)

	# set up output
	out <- list()
	out[["param"]] <- param

	# check Scans
	scan.include <- getScanIndex(data = seqData, scan.include = nullModObj$scanID)	
	# set up a filter
	seqSetFilter(seqData, sample.sel = scan.include$index, verbose=FALSE)
	
	# logical of samples to use for Allele Frequency calculations
	if(is.null(AF.sample)){
		AF.sample.use <- rep(TRUE, scan.include$n)
	}else{
		AF.sample.use <- scan.include$value %in% AF.sample
	}
	nsample <- list()
	nsample[["analysis"]] <- scan.include$n
	nsample[["AF"]] <- sum(AF.sample.use)
	out[["nsample"]] <- nsample
	
	# check that all variants in aggVarList are in seqData
	variant.include <- as.vector(unlist(sapply(aggVarList, function(x){ x$variant.id }), use.names=FALSE))
	variant.include <- getVariantInclude(data = seqData, variant.include = variant.include, chromosome = NULL)
	

	# set up main results matrix
	nv <- c("n.site", "n.sample.alt")
        nv <- .outputColumns(nv, AF.range = AF.range, weight.beta = weight.beta, weight.user = weight.user,
                             test = test, burden.test = burden.test, rho = rho, verbose=verbose)

	# determine the number of variant blocks
	nblocks <- length(aggVarList)
	if(verbose) message("Running analysis with ", scan.include$n, " Samples and ", variant.include$n, " Variants in ", nblocks, " Aggregate Units")

	resMain <- matrix(NA, nrow=nblocks, ncol=length(nv), dimnames=list(NULL,nv))

	# set up results for variants
	nvm <- c("variantID", "allele", "chr", "pos", "n.obs", "freq", "weight")
	variantInfo <- list()


	# get residuals and sqrt of Projection matrix
	proj <- .calculateProjection(nullModObj = nullModObj, test = test, burden.test = burden.test)

	# keep track of time for rate reporting
	startTime <- Sys.time()


	# loop through aggVarList
	for(b in 1:nblocks){		

		# extract variant info for the block
		aggVarBlock <- aggVarList[[b]]

                # subset aggVarBlock to requested/polymorphic variants
                aggVarBlock <- aggVarBlock[aggVarBlock$variant.id %in% variant.include$value,]
                
		# set a filter for variants in this group
                variant.index <- variant.include$index[variant.include$value %in% aggVarBlock$variant.id]
		seqSetFilter(seqData, variant.sel = variant.index, verbose = FALSE)

		# get list of allele indices to test for each variant
		alleleIndex <- do.call(list, by(aggVarBlock, aggVarBlock$variant.id, function(x){ x$allele.index }, simplify = FALSE))
		# number of alleles per variant
		nAllele <- sapply(alleleIndex, length)

		# matrix to store variant information for this group
		variantRes <- matrix(NA, nrow = nrow(aggVarBlock), ncol = length(nvm), dimnames = list(NULL,nvm))
                if (nrow(aggVarBlock) == 0) {
                        resMain[b,"n.site"] <- 0
			variantInfo[[b]] <- as.data.frame(variantRes)
			next
                }
		
		variantRes[,"variantID"] <- aggVarBlock$variant.id
		variantRes[,"allele"] <- aggVarBlock$allele.index
		# read in chr, pos (repeat the appropriate number of times for each variant)
		chromChar <- rep(seqGetData(seqData, "chromosome"), nAllele)
		#variantRes[,"chr"] <- as.numeric(chromChar)
		variantRes[,"pos"] <- rep(seqGetData(seqData, "position"), nAllele)
                
		# get genotypes for the block; rows are samples, columns are variant-alleles
		geno <- do.call(cbind, alleleDosage(seqData, n = alleleIndex))
		
		# calculate number of missing/observed values by variant-allele
		n.miss <- colSums(is.na(geno))
		variantRes[,"n.obs"] <- nrow(geno) - n.miss

		# calculate allele frequencies
                sex <- sampleData(seqData)$sex
		freq <- alleleFreq(geno = geno, chromChar = chromChar, sex = sex, sample.use = AF.sample.use)
		variantRes[,"freq"] <- freq

		# index of those inside the allele freq threshold
		include <- (freq !=0 & freq !=1 & freq >= AF.range[1] & freq <= AF.range[2])
		variantRes[!include,"weight"] <- 0		
		
		# number of variant sites
		n.site <- length(unique(variantRes[include,"variantID"]))
		resMain[b,"n.site"] <- n.site

		if(n.site == 0){	
			variantInfo[[b]] <- as.data.frame(variantRes)
                        variantInfo[[b]][,"chr"] <- chromChar
			next

		}else{
			# subset
			geno <- geno[, include, drop = FALSE]
			freq <- freq[include]

			# calculate number of samples with observed alternate alleles > 0
			resMain[b,"n.sample.alt"] <- sum(rowSums(geno, na.rm=TRUE) > 0)

			# variant weights
			if(is.null(weight.user)){
				#freq <- ifelse(freq < 0.5, freq, 1-freq) - should MAF be used even if alternate allele is not minor?
				# Beta weights
				weights <- dbeta(freq, weight.beta[1], weight.beta[2])
			}else{
				# user supplied weights - read in from variantData, repeat to match nAllele, subset to those included
				weights <- rep(pData(variantData(seqData))[,weight.user], nAllele)[include]
			}
			variantRes[include,"weight"] <- weights
			variantInfo[[b]] <- as.data.frame(variantRes)
                        variantInfo[[b]][,"chr"] <- chromChar
		}

		# mean impute missing genotype values
		if(sum(n.miss) > 0){
			miss.idx <- which(is.na(geno))
			miss.snp.idx <- ceiling(miss.idx/nrow(geno))
			geno[miss.idx] <- 2*freq[miss.snp.idx]
		}

		# perform test
		if(test == "Burden"){
			# multiply each variant by weights and collapse
			burden <- colSums(t(geno)*weights)
			testout <- .runBurdenTest(burden = burden, projObj = proj, burden.test = burden.test)

		}else if(test == "SKAT"){
			U <- as.vector(crossprod(geno, proj$resid))
			geno <- crossprod(proj$Mt, geno)
			if(length(rho) == 1){
				testout <- .runSKATTest(scores = U, geno.adj = geno, weights = weights, rho = rho, pval.method = pval.method, optimal = FALSE)
			}else{
				testout <- .runSKATTest(scores = U, geno.adj = geno, weights = weights, rho = rho, pval.method = pval.method, optimal = TRUE)
			}
		}

		# update main results
		for(val in names(testout)){ resMain[b,val] <- testout[val] }

		# report time				
		if(verbose & b %% 100 == 0){
			endTime <- Sys.time()
			rate <- format(endTime - startTime, digits=4)
			message(paste("...Aggregate Block", b, "of", nblocks, "Completed -", rate, "Elapsed"))
		}		
	}

	# add results to main output
	resMain <- as.data.frame(resMain)
	rownames(resMain) <- names(aggVarList)	
	out[["results"]] <- resMain
	# add variantInfo to main output
	names(variantInfo) <- names(aggVarList)
	out[["variantInfo"]] <- variantInfo

	# return to original filter
	seqSetFilter(seqData, sample.sel = seqFilt.original$sample.sel, variant.sel = seqFilt.original$variant.sel, verbose = FALSE)
	
	return(out)	
}



assocTestSeqWindow <- function(seqData,
                               nullModObj,
                               variant.include = NULL,
                               chromosome = NULL,
                               window.size = 50,
                               window.shift = 20,
                               AF.sample = NULL,
                               AF.range = c(0,1),
                               weight.beta = c(0.5, 0.5),
                               weight.user = NULL,
                               test = "Burden",
                               burden.test = "Score",
                               rho = 0,
                               pval.method = "kuonen",
                               verbose = TRUE){

	# save the filter
	seqFilt.original <- seqGetFilter(seqData)
        # reset so indexing works
        seqResetFilter(seqData, verbose=FALSE)

	# check the parameters
	param <- .paramChecks(seqData = seqData, AF.range = AF.range, weight.beta = weight.beta, weight.user = weight.user, 
                              test = test, burden.test = burden.test, family = nullModObj$family$family, 
                              mixedmodel = nullModObj$family$mixedmodel, rho = rho, pval.method = pval.method)

	# set up output
	out <- list()
	out[["param"]] <- param

	# checks on windows
	if(window.size <= 0 | window.shift <= 0){
		stop("window.size and window.shift must be positive")
	}
	if(window.size < window.shift){
		stop("window.shift can not be larger than window.size")
	}
	window.param <- list()
	window.param[["size"]] <- window.size
	window.param[["shift"]] <- window.shift
	out[["window"]] <- window.param
	
	# check Scans
	scan.include <- getScanIndex(data = seqData, scan.include = nullModObj$scanID)	
	# set up a filter
	seqSetFilter(seqData, sample.sel = scan.include$index, verbose=FALSE)	

	# logical of samples to use for Allele Frequency calculations
	if(is.null(AF.sample)){
		AF.sample.use <- rep(TRUE, scan.include$n)
	}else{
		AF.sample.use <- scan.include$value %in% AF.sample
	}
	nsample <- list()
	nsample[["analysis"]] <- scan.include$n
	nsample[["AF"]] <- sum(AF.sample.use)
	out[["nsample"]] <- nsample

	# check Variants and set filter - filters monomorphic reference alleles
	variant.include <- getVariantInclude(data = seqData, variant.include = variant.include, chromosome = chromosome)
	# set a filter
        seqSetFilter(seqData, variant.sel = variant.include$index, verbose = FALSE)
        # get chromosome
        variant.chr <- seqGetData(seqData, "chromosome")


	# set up main results matrix
	nv <- c("chr", "window.start", "window.stop", "n.site", "dup")
        nv <- .outputColumns(nv, AF.range = AF.range, weight.beta = weight.beta, weight.user = weight.user,
                             test = test, burden.test = burden.test, rho = rho, verbose=verbose)
	resMain <- matrix(NA, nrow = 0, ncol = length(nv), dimnames = list(NULL,nv))

	# set up results for variants
	nvm <- c("variantID", "allele", "chr", "pos", "n.obs", "freq", "weight")
	variantInfo <- matrix(NA, nrow=0, ncol=length(nvm), dimnames = list(NULL,nvm))


	# get residuals and sqrt of Projection matrix
	proj <- .calculateProjection(nullModObj = nullModObj, test = test, burden.test = burden.test)
	
	# loop through chromosomes
	for(chr in unique(variant.chr)){
		# keep track of time for rate reporting
		startTime <- Sys.time()

		# variant index for this chromosome
		variantID <- variant.include$index[variant.chr == chr]
		# set up a filter for the chromosome
		seqSetFilter(seqData, variant.sel = variantID, verbose = FALSE)
		# get position
		variantPos <- seqGetData(seqData, "position")

		# idenfity windows
		window.start <- seq(from = (ceiling(variantPos[1]/(window.shift*1000))-1)*window.shift*1000 + 1,
							to = tail(variantPos,1) - window.size*1000 + window.shift*1000,
							by = window.shift*1000)
		window.stop <- window.start + window.size*1000 - 1
		# determine the number of variant windows
		nblocks <- length(window.start)

		# results for this chromosome
		resChr <- matrix(NA, nrow = nblocks, ncol = length(nv), dimnames = list(NULL,nv))
		resChr[,"chr"] <- .chrToInt(chr)

		if(verbose) message("Running analysis on Chromosome ", chr, " with ", scan.include$n, " Samples and ", length(variantID), " Variants in ", nblocks, " Sliding Windows of size ", window.size, "kb")


		# set some initial parameters
		count <- 0
		testID <- NULL
		prevWindowID <- NULL
		geno <- NULL
		weights <- NULL
		if(test == "Burden"){
			burden <- numeric(scan.include$n)
		}else if(test == "SKAT"){
			U <- NULL
		}

		# slide through the windows
		for(b in 1:nblocks){
			
			# IDs of variants in the window
			windowID <- variantID[variantPos >= window.start[b] & variantPos <= window.stop[b]]
			
			# check if there are any variants
			if(length(windowID) > 0){
				count <- count + 1 
				resChr[count, "window.start"] <- window.start[b]
				resChr[count, "window.stop"] <- window.stop[b]

				# check if there are variants to drop
				var.drop <- !(testID %in% windowID)
				if(any(var.drop)){
					# update
					if(test == "Burden"){
						# remove burden of variants to drop
						burden.drop <- colSums(t(geno[,var.drop])*weights[var.drop])
						burden <- burden - burden.drop
					}else if(test == "SKAT"){
						# remove scores of variants to drop
						U <- U[!var.drop]
					}

					# remove genotypes of variants to drop
					geno <- geno[,!var.drop, drop=FALSE]
					# remove weights of variants to drop
					weights <- weights[!var.drop]
					# remove variantIDs of variants to drop
					testID <- testID[!var.drop]
				}

				# check if there are any variants to add
				var.add <- windowID[!(windowID %in% prevWindowID)]
				if(length(var.add) > 0){
					# set a filter for variants to add
					seqSetFilter(seqData, variant.sel = var.add, verbose = FALSE)

					# number of alleles at each variant
					nAllele <- seqNumAllele(seqData) - 1
					# get list of allele indices to test for each variant
					alleleIndex <- list(); for(i in 1:length(nAllele)){ alleleIndex[[i]] <- 1:nAllele[i] }

					# matrix to store variant information					
					variantRes <- matrix(NA, nrow=sum(nAllele), ncol=length(nvm), dimnames = list(NULL,nvm))
					# repeat each value the appropriate number of times for each variant
                                        var.id <- variant.include$value[variant.include$index %in% var.add]
					variantRes[,"variantID"] <- rep(var.id, nAllele)
					variantRes[,"allele"] <- unlist(alleleIndex)
					chromChar <- rep(chr, sum(nAllele))
					variantRes[,"chr"] <- .chrToInt(chromChar)
					variantRes[,"pos"] <- rep(variantPos[variantID %in% var.add], nAllele)

					# get genotypes; rows are samples, columns are variant-alleles
					geno.add <- do.call(cbind, alleleDosage(seqData, n = alleleIndex))

					# calculate number of missing/observed values by variant-allele
					n.miss <- colSums(is.na(geno.add))
					variantRes[,"n.obs"] <- nrow(geno.add) - n.miss

					# calculate allele frequencies
                                        sex <- sampleData(seqData)$sex
					freq <- alleleFreq(geno = geno.add, chromChar = chromChar, sex = sex, sample.use = AF.sample.use)
					variantRes[,"freq"] <- freq

					# index of those inside the allele freq threshold
					include <- (freq !=0 & freq !=1 & freq >= AF.range[1] & freq <= AF.range[2])
					variantRes[!include,"weight"] <- 0

					# subset
					geno.add <- geno.add[,include, drop=FALSE]
					freq <- freq[include]

					# variant weights
					if(is.null(weight.user)){
						#freq <- ifelse(freq < 0.5, freq, 1-freq) - should MAF be used even if alternate allele is not minor?
						# Beta weights
						weights.add <- dbeta(freq, weight.beta[1], weight.beta[2])
					}else{
						# user supplied weights - read in from variantData, repeat to match nAllele, subset to those included
						weights.add <- rep(pData(variantData(seqData))[,weight.user], nAllele)[include]
					}
					variantRes[include,"weight"] <- weights.add

					# mean impute missing genotype values
					if(sum(n.miss) > 0){
						miss.idx <- which(is.na(geno.add))
						miss.snp.idx <- ceiling(miss.idx/nrow(geno.add))
						geno.add[miss.idx] <- 2*freq[miss.snp.idx]
					}

					# update 
					variantInfo <- rbind(variantInfo, variantRes)

					if(test == "Burden"){
						# add burden of variants to add
						burden <- burden + colSums(t(geno.add)*weights.add)
						# add genotypes of variants to add
						geno <- cbind(geno, geno.add)
					}else if(test == "SKAT"){
						# add scores of variants to add
						U <- c(U, as.vector(crossprod(geno.add, proj$resid)))
						# add genotypes of variants to add
						geno <- cbind(geno, crossprod(proj$Mt, geno.add))
					}
					# add weights of variants to add
					weights <- c(weights, weights.add)

					testID <- c(testID, variant.include$index[variant.include$value %in% variantRes[include,"variantID"]])
				}

				# number of variant sites
				n.site <- length(unique(testID))
				resChr[count,"n.site"] <- n.site

				if(n.site > 0){
					if(any(var.drop) | length(var.add) > 0){
						resChr[count,"dup"] <- 0
						# perform test
						if(test == "Burden"){
							testout <- .runBurdenTest(burden = burden, projObj = proj, burden.test = burden.test)

						}else if(test == "SKAT"){
							if(length(rho) == 1){
								testout <- .runSKATTest(scores = U, geno.adj = geno, weights = weights, rho = rho, pval.method = pval.method, optimal = FALSE)
							}else{
								testout <- .runSKATTest(scores = U, geno.adj = geno, weights = weights, rho = rho, pval.method = pval.method, optimal = TRUE)
							}
						}

					}else{
						resChr[count,"dup"] <- 1
					}
					# update main results
					for(val in names(testout)){ resChr[count,val] <- testout[val] }
				}

			# no variants in the window
			}else{ 
				# drop a line from the results
				resChr <- resChr[-(count+1),]
			}

			# record this window's variantIDs
			prevWindowID <- windowID

			# report time				
			if(verbose & b %% 100 == 0){
				endTime <- Sys.time()
				rate <- format(endTime - startTime, digits=4)
				message(paste("...Window", b, "of", nblocks, "Completed -", rate, "Elapsed"))
			}

		} # for windows
		# add results to main set
		resMain <- rbind(resMain, resChr)	

	} # for chromosomes

	# add results to main output
	out[["results"]] <- as.data.frame(resMain)
        out[["results"]][,"chr"] <- .intToChr(out[["results"]][,"chr"])
	# add variantInfo to main output	
	out[["variantInfo"]] <- as.data.frame(variantInfo)
        out[["variantInfo"]][,"chr"] <- .intToChr(out[["variantInfo"]][,"chr"])

	# return to original filter
	seqSetFilter(seqData, sample.sel = seqFilt.original$sample.sel, variant.sel = seqFilt.original$variant.sel, verbose = FALSE)
	
	return(out)
}









.paramChecks <- function(seqData, AF.range, weight.beta, weight.user, test, burden.test, family, mixedmodel, rho, pval.method){
	# list of parameters
	param <- list()
	
	# check AF.range
	if(AF.range[1] < 0 | AF.range[2] > 1){ stop("AF.range must be in [0,1]")}
	param[["AF.range"]] <- AF.range

	# check weights
	if(is.null(weight.user)){
		if(length(weight.beta) != 2){ stop("weight.beta must be a vector of length 2 specifying the Beta parameters for the weight function")}
		param[["weight.beta"]] <- weight.beta
		param[["weight.user"]] <- FALSE
	}else{
		param[["weight.beta"]] <- FALSE
		if(!(weight.user %in% varLabels(variantData(seqData)))){ stop("The variable specified for weight.user must be in the variantData slot of seqData") }
		param[["weight.user"]] <- weight.user
	}

	param[["family"]] <- family
	param[["mixedmodel"]] <- mixedmodel
	
	# check that test is valid
	if(!is.element(test,c("Burden","SKAT"))){ stop("test must be Burden or SKAT") }
	param[["test"]] <- test
	if(test == "Burden"){
		# check burden.test type
		if(family == "gaussian"){
			if(!is.element(burden.test,c("Score","Wald"))){ stop("burden.test must be one of Score or Wald when family is gaussian")}
		}
		if(!mixedmodel & family == "binomial"){
			if(!is.element(burden.test,c("Score","Wald","Firth"))){ stop("burden.test must be one of Score, Wald, or Firth for a glm with binomial family")}
			if(burden.test == "Firth") requireNamespace("logistf")
		}
		if(mixedmodel & family == "binomial"){
			if(!is.element(burden.test,c("Score"))){ stop("burden.test must be Score for a mixed model with binomial family")}
		}
		param[["burden.test"]] <- burden.test		
	}	
	if(test == "SKAT"){
		# check rho value
		if(any(rho < 0) | any(rho > 1)){ stop("rho must take values in [0,1]") }
		if(length(rho) > 1){ param[["test"]] <- "SKAT-O" }
		# check pval.method
		if(!(pval.method %in% c("kuonen","davies","liu"))){ stop("pval.method must be one of 'kuonen', 'davies', or 'liu'")}
		if(pval.method == "kuonen") requireNamespace("survey")
		if(pval.method == "davies" | pval.method == "liu") requireNamespace("CompQuadForm")
		param[["rho"]] <- rho
		param[["pval.method"]] <- pval.method
	}
	
	return(param)	
}



.outputColumns <- function(nv, AF.range, weight.beta, weight.user, test, burden.test, rho, verbose){
	if(test == "Burden"){
		nv <- append(nv, "burden.skew")
		if(burden.test == "Score"){		
			nv <- append(nv, c("Score", "Var", "Score.stat", "Score.pval"))
		}else if(burden.test == "Wald"){
			nv <- append(nv, c("Est", "SE", "Wald.stat", "Wald.pval"))
		}else if(burden.test == "Firth"){
			nv <- append(nv, c("Est", "SE", "Firth.stat", "Firth.pval"))
		}
		if(verbose){
			if(is.null(weight.user)){
				message("Performing ", burden.test, " Burden Tests for Variants with AF in [", AF.range[1], ",", AF.range[2], "] using weights from a Beta(", weight.beta[1], ",", weight.beta[2], ") distribution")
			}else{
				message("Performing ", burden.test, " Burden Tests for Variants with AF in [", AF.range[1], ",", AF.range[2], "] using weights specified by ", weight.user, " in the variantData slot of seqData")			
			}
		}	
		
	}else if(test == "SKAT"){
		nv <- append(nv, c(paste("Q",rho,sep="_"), paste("pval",rho,sep="_"), paste("err",rho,sep="_")))
		if(length(rho) == 1){			
			if(verbose){
				if(is.null(weight.user)){
					message("Performing SKAT Tests for Variants with AF in [", AF.range[1], ",", AF.range[2], "] using weights from a Beta(", weight.beta[1], ",", weight.beta[2], ") distribution and rho = ", rho)
				}else{
					message("Performing SKAT Tests for Variants with AF in [", AF.range[1], ",", AF.range[2], "] using weights specified by ", weight.user, " in the variantData slot of seqData and rho = ", rho)
				}
			}
		}else{
			nv <- append(nv, c("min.pval", "opt.rho", "pval_SKATO"))
			if(verbose){
				if(is.null(weight.user)){
					message("Performing SKAT-O Tests for Variants with AF in [", AF.range[1], ",", AF.range[2], "] using weights from a Beta(", weight.beta[1], ",", weight.beta[2], ") distribution and rho = (", paste(rho, collapse=", "), ")")
				}else{
					message("Performing SKAT-O Tests for Variants with AF in [", AF.range[1], ",", AF.range[2], "] using weights specified by ", weight.user, " in the variantData slot of seqData and rho = (", paste(rho, collapse=", "), ")")
				}
			}
		}		
	}

        return(nv)
}



.calculateProjection <- function(nullModObj, test, burden.test){
	# covariate matrix
	W <- nullModObj$model.matrix
	# outcome
	Y <- nullModObj$workingY

	if(nullModObj$family$mixedmodel){
		C <- nullModObj$cholSigmaInv
		CW <- crossprod(C,W)
		# Projection matrix P = Mt %*% M
		Mt <- C - tcrossprod(tcrossprod(C,tcrossprod(chol2inv(chol(crossprod(CW))),CW)),CW); rm(C); rm(CW)
		# residuals
		resid <- as.vector(Mt %*% crossprod(Mt,Y))
		return(list(W = W, Y = Y, Mt = Mt, resid = resid, family = nullModObj$family$family))

	}else if(nullModObj$family$family == "gaussian"){
		# Projection matrix P = Mt %*% M
		Mt <- (diag(nrow(W)) - tcrossprod(tcrossprod(W,chol2inv(chol(crossprod(W)))),W))/nullModObj$sigma
		# residuals
		resid <- nullModObj$resid.response/(nullModObj$sigma^2)
		return(list(W = W, Y = Y, Mt = Mt, resid = resid, family = nullModObj$family$family))

	}else if(nullModObj$family$family == "binomial" & !(test == "Burden" & (burden.test == "Wald" | burden.test == "Firth")) ){
		sigma <- nullModObj$sigma
		C <- diag(sigma)
		CW <- W*sigma
		# Projection matrix P = Mt %*% M
		Mt <- C - tcrossprod(tcrossprod(C,tcrossprod(chol2inv(chol(crossprod(CW))),CW)),CW); rm(C); rm(CW)
		# residuals
		resid <- nullModObj$resid.response
		return(list(W = W, Y = Y, Mt = Mt, resid = resid, family = nullModObj$family$family))

	}else{
		return(list(W = W, Y = Y, family = nullModObj$family$family))
	}
}



.runBurdenTest <- function(burden, projObj, burden.test){
	# calculate skewness of the burden
	burden.skew <- .skewness(burden)

	if(burden.test == "Score"){
		Xtilde <- crossprod(projObj$Mt,burden)
		XtX <- sum(Xtilde^2)
		XtY <- sum(burden*projObj$resid)
		stat <- (XtY)^2/XtX
		return(c(burden.skew = burden.skew, Score = XtY, Var = XtX, Score.stat = stat, Score.pval = pchisq(stat, df=1, lower.tail=FALSE)))

	}else if(burden.test == "Wald"){
		if(projObj$family == "gaussian"){
			Xtilde <- crossprod(projObj$Mt,burden)
			XtX <- sum(Xtilde^2)
			XtY <- sum(burden*projObj$resid)
			beta <- XtY/XtX
			RSS <- (sum(projObj$Y*projObj$resid) - XtY*beta)/(length(projObj$Y) - ncol(projObj$W) - 1)
			Vbeta <- RSS/XtX
			stat <- beta^2/Vbeta
			return(c(burden.skew = burden.skew, Est = beta, SE = sqrt(Vbeta), Wald.stat = stat, Wald.pval = pchisq(stat, df=1, lower.tail=FALSE)))
			# beta and SE match what is given by lm() when not a mixed model; this is a Z test, lm() performs a t test

		}else if(projObj$family == "binomial"){
			vals <- c(burden.skew, summary(glm(projObj$Y ~ projObj$W + burden - 1, family = "binomial"))$coef["burden",])
			names(vals) <- c("burden.skew", "Est", "SE", "Wald.stat", "Wald.pval")
			vals["Wald.stat"] <- vals["Wald.stat"]^2		
			return(vals)
		}

	}else if(burden.test == "Firth"){
		modfit <- logistf::logistf(projObj$Y ~ projObj$W + burden - 1)
		idx <- which(names(modfit$coef) == "burden")
		pval <- modfit$prob["burden"]
		vals <- c(burden.skew, modfit$coef[idx], sqrt(modfit$var[idx,idx]), qchisq(pval, df=1, lower.tail = FALSE), pval)
		names(vals) <- c("burden.skew", "Est", "SE", "Firth.stat", "Firth.pval")
		return(vals)
	}
}

.skewness <- function(x){
	mu3 <- mean( (x - mean(x))^3 )
	sigma3 <- var(x)^(3/2)
	mu3/sigma3
}



.runSKATTest <- function(scores, geno.adj, weights, rho, pval.method, optimal){
	# covariance of scores
	V <- crossprod(geno.adj)
	
	# vector to hold output
	out <- numeric(length = 3*length(rho))
	names(out)  <- c(paste("Q",rho,sep="_"), paste("pval",rho,sep="_"), paste("err",rho,sep="_"))

	# for SKAT-O need to save lambdas
	if(optimal){ lambdas <- vector("list",length(rho)) }

	for(i in 1:length(rho)){
		if(rho[i] == 0){
			# Variance Component Test
			Q <- sum((weights*scores)^2) # sum[(w*scores)^2]  # for some reason SKAT_emmaX divides this by 2
			distMat <- weights*t(weights*V)  # (weights) V (weights) = (weights) X' P X (weights)

		}else if(rho[i] == 1){
			# Burden Test
			Q <- sum(weights*scores)^2 # (sum[w*scores])^2  # for some reason SKAT_emmaX divides this by 2
			distMat <- crossprod(weights,crossprod(V, weights)) # weights^T V weights

		}else if(rho[i] > 0 & rho[i] < 1){
			rhoMat <- matrix(rho[i], nrow=length(scores), ncol=length(scores)); diag(rhoMat) <- 1
			cholRhoMat <- t(chol(rhoMat, pivot=TRUE))
			Q <- crossprod(crossprod(weights*cholRhoMat,scores)) # scores' (weights) (rhoMat) (weights) scores
			distMat <- crossprod(cholRhoMat, crossprod(weights*t(weights*V), cholRhoMat)) # (cholRhoMat) (weights) X' P X (weights) (cholRhoMat)
		}

		# lambda for p value calculation
		lambda <- eigen(distMat, only.values = TRUE, symmetric=TRUE)$values
		lambda <- lambda[lambda > 0]
		if(optimal){ lambdas[[i]] <- lambda }

		# p value calculation
		if(length(scores) == 1){
			pval <- pchisq(Q/distMat, df=1, lower.tail=FALSE)
			err <- 0

		}else{
			if(pval.method == "kuonen"){
				pval <- survey:::saddle(x = Q, lambda = lambda)
				err <- ifelse(is.na(pval), 1, 0)

			}else if(pval.method == "davies"){
				tmp <- CompQuadForm::davies(q = Q, lambda = lambda, acc = 1e-06)
				pval <- tmp$Qq
				err <- tmp$ifault

			}else if(pval.method == "liu"){
				pval <- CompQuadForm::liu(q = Q, lambda = lambda)
				err <- 0
			}

			if(err > 0){
				pval <- CompQuadForm::liu(q = Q, lambda = lambda)
			}
		}

		# update results
		out[c(paste("Q",rho[i],sep="_"), paste("pval",rho[i],sep="_"), paste("err",rho[i],sep="_"))] <- c(Q, pval, err)
	}

	# SKAT-O
	if(optimal){
		# vector for output
		out2 <- rep(NA, 3)
		names(out2) <- c("min.pval", "opt.rho", "pval_SKATO")

		if(length(scores) == 1){
			# pvalue is the same for all rhos
			pval <- out[grep("pval", names(out))][1]
			out2["min.pval"] <- pval
			out2["pval_SKATO"] <- pval

		}else{
			# find the minimum p-value
			pval <- out[grep("pval", names(out))]
			minp <- min(pval)
			out2["min.pval"] <- minp
			out2["opt.rho"] <- rho[which.min(pval)]

			# get qmin(rho); i.e. the (1-minp)th percentile of dist of each Q
			qmin <- rep(NA, length(rho))
			for(i in 1:length(rho)){
				qmin[i] <- skatO_qchisqsum(minp, lambdas[[i]])
			}

			# calculate other terms
			Z <- t(t(geno.adj)*weights)
			zbar <- rowMeans(Z)
			zbarTzbar <- sum(zbar^2)
			M <- tcrossprod(zbar)/zbarTzbar
			ZtImMZ <- crossprod(Z, crossprod(diag(nrow(M)) - M, Z))
			lambda.k <- eigen(ZtImMZ, symmetric = TRUE, only.values = TRUE)
			lambda.k <- lambda.k$values[lambda.k$values > 0]
			mua <- sum(lambda.k)
			sum.lambda.sq <- sum(lambda.k^2)
			sig2a <- 2*sum.lambda.sq
			trMatrix <- crossprod(crossprod(Z,crossprod(M,Z)),ZtImMZ)
			sig2xi <- 4*sum(diag(trMatrix))
			kera <- sum(lambda.k^4)/sum.lambda.sq^2 * 12
			ldf <- 12/kera

			# calculate tau(rho)
			tau <- ncol(Z)^2*rho*zbarTzbar + (1-rho)*sum(crossprod(zbar, Z)^2)/zbarTzbar

			# find min{(qmin(rho)-rho*chisq_1)/(1-rho)} with integration							
			otherParams <- c(mu = mua, degf = ldf, varia = sig2a+sig2xi)
			# integrate
			re <- tryCatch({
                            integrate(integrateFxn, lower = 0, upper = 40, subdivisions = 2000, qmin = qmin, otherParams = otherParams, tau = tau, rho = rho, abs.tol = 10^-25)
                        }, error=function(e) NA)
			out2["pval_SKATO"] <- 1-re[[1]]
		}
		
		# update results
		out <- append(out, out2)
	}

	# return results
	return(out)	
}


# function to calculate q_min value
# basically a qchisqsum() function that takes the quantile/percentile and the lambda values
# matches the first 2 moments and the kurtosis
# based upon liu et al (2009) paper
skatO_qchisqsum <- function(p, lambdas){
	mu <- sum(lambdas)
	sum.lambda.sq <- sum(lambdas^2)

  	s1 <- sum(lambdas^3)/(sum.lambda.sq^(3/2))
  	s2 <- sum(lambdas^4)/(sum.lambda.sq^2)
  	if(s1^2 > s2){
    	a <- 1/(s1-sqrt(s1^2-s2))
    	d <- s1*a^3 - a^2
    	l <- a^2 - 2*d
  	}else{ # s1^2 <= s2
            l <- 1/s2 # in liu et al, this is l=1/s1^2; matches kurtosis instead of skewness to improve tail prob estimates
        }  	
  
  	qmin <- qchisq(1-p, df=l)
  	pval <- (qmin - l)/sqrt(2*l) * sqrt(2*sum.lambda.sq) + mu
  
  	return(pval)
}


## function to integrate; the first term of the optimal integrand
# it's a non-central sum of weighted chi-squares
integrateFxn <- function(x, qmin, otherParams, tau, rho){
    n.r <- length(rho)
    n.x <- length(x)
    
    t1 <- tau %x% t(x)
    tmp <- (qmin - t1)/(1-rho)
    minval <- apply(tmp,2,min)

    degf <- otherParams["degf"]
    mu <- otherParams["mu"]
    varia <- otherParams["varia"]

    temp.q<-(minval - mu)/sqrt(varia)*sqrt(2*degf) + degf

    re<-pchisq(temp.q ,df=degf) * dchisq(x,df=1)

    return(re)
}


.chrToInt <- function(chr){
    unname(setNames(1:24, c(1:22, "X", "Y"))[chr])
}

.intToChr <- function(chr){
    unname(setNames(c(1:22, "X", "Y"), 1:24)[as.character(chr)])
}
