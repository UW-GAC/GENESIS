fitNullReg <- function(	scanData,
						outcome,
						covars = NULL,
						scan.include = NULL,
						family = gaussian,
						verbose = TRUE){

	# make sure scanData is a data.frame
	if(class(scanData) == "ScanAnnotationDataFrame"){
		scanData <- pData(scanData)
	}
	if(class(scanData) == "AnnotatedDataFrame"){
		scanData <- pData(scanData)
		names(scanData)[names(scanData) == "sample.id"] <- "scanID"
	}
	if(class(scanData) != "data.frame"){
		stop("scanData should either be a data.frame, an AnnotatedDataFrame, or a ScanAnnotationDataFrame from GWASTools")
	}

	# check family
    if(is.character(family)){
        family <- get(family)
    }
    if(is.function(family)){
        family <- family()
    }
    if(is.null(family$family)){
        stop("'family' not recognized")
    }

    # create design matrices and outcome vector
    if(verbose) message("Reading in Phenotype and Covariate Data...")
    dat <- createDesignMatrix(scanData = scanData, outcome = outcome, covars = covars, scan.include = scan.include)

    # update scan.include
    scan.include <- getScanIndex(data = scanData, scan.include = rownames(dat$W))    
    if(verbose) message("Fitting Model with ", scan.include$n, " Samples")


    # fit the null model
    if(family$family == "gaussian"){
    	fit0 <- lm(dat$Y ~ -1 + dat$W)
    	# sqrt(RSS)
    	sigma <- summary(fit0)$sigma
    	# residuals
		resid <- residuals(fit0, type="response")

	}else if(family$family == "binomial"){
		fit0 <- glm(dat$Y ~ -1 + dat$W, family=family)
		# variance function
		sigma <- sqrt(family$variance(fit0$fitted))  # mu(1-mu) for binomial
		# residuals
		resid <- residuals(fit0, type="response")
	}

	fixef <- as.data.frame(summary(fit0)$coef)
	varNames <- gsub(pattern = "dat[$]W", replacement = "", x = rownames(fixef))
	rownames(fixef) <- varNames

	betaCov <- vcov(fit0)
	dimnames(betaCov) <- list(varNames, varNames)

	# check for aliased (removed) coefficients and remove from dat$W
	aliased <- summary(fit0)$aliased
	names(aliased) <- colnames(dat$W)
	dat$W <- dat$W[,!aliased, drop = FALSE]

	family$mixedmodel <- FALSE

	return(list(fixef = fixef,
				betaCov = betaCov,
				fitted.values = fit0$fitted.values,
				resid.response = resid,
				logLik = as.numeric(logLik(fit0)),
				AIC = AIC(fit0),
				workingY = as.vector(dat$Y),
				model.matrix = dat$W,
				aliased = aliased,
				sigma = sigma,
				scanID = scan.include$value,
				family = family))

}
