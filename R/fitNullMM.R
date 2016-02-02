fitNullMM <- function(	scanData,						
						outcome, 
						covars = NULL,
						covMatList,
						scan.include = NULL,						 
						family = gaussian, 
						group.var = NULL,
						start = NULL,
						AIREML.tol = 1e-6,
						maxIter = 100,
						dropZeros = TRUE,
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
		stop("scanData should either be a data.frame, and AnnotatedDataFrame, or a ScanAnnotationDataFrame from GWASTools")
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

    # determine unique groups and create index
    if(!is.null(group.var)){
        if(family$family != "gaussian"){
            stop("heterogeneous group residual variances can only be used when family = gaussian")
        }
        if(!(group.var %in% names(scanData))){
        	stop("group.var must be in scanData")
        }
        # read in group variable and subset to included samples
        group <- scanData[scan.include$index, group.var]
        # determine unique groups
        group.names <- as.character(unique(group))
        # number of groups
        g <- length(group.names)

        # index of samples in each group
        group.idx <- list()
        for(i in 1:g){
            group.idx[[i]] <- which(group == group.names[i])
        }

    }else{
    	if(family$family == "gaussian"){
    		group.names <- "E"
    		g <- 1
    	}else{
    		g <- 0
    	}
    }
    
    # if covMatList is a matrix, convert to a list
    if(class(covMatList) == "matrix"){
        covMatList <- list(A = covMatList)
    }

    # number of covariance structure matrices
    m <- length(covMatList)

    # if covMatList doesn't have names, assign them
    if(is.null(names(covMatList))){
        names(covMatList) <- paste("A",1:m,sep="")
    }
    
    # check for starting values
    if(!is.null(start)){
        if(family$family == "gaussian"){
            if(length(start) != (m+g)){
                stop("length of start must equal the number of matrices in covMatList + number of groups for gaussian")
            }
        }else{
            if(length(start) != m){
                stop("length of start must equal the number of matrices in covMatList")
            }
        }
    }
    
    # subest covariance structure matrices
    for(i in 1:m){
        if(!all(scan.include$value %in% colnames(covMatList[[i]]))){
            stop(paste("All of the included Samples must be in matrix", i, "of covMatList"))
        }
        # subset matrix
        keepMat <- colnames(covMatList[[i]]) %in% scan.include$value
        covMatList[[i]] <- covMatList[[i]][keepMat,keepMat]
        
        # check that names match
        if(!all(colnames(covMatList[[i]]) == scan.include$value)){
            stop("Column and Row names of matrix", i, "of covMatList must match the scanID in scanData")
        }
    }
    

    if(verbose) message("Computing Variance Component Estimates using AIREML Procedure...")
        
    if(family$family == "gaussian"){
        # estimate variance components        
        if(verbose) message(paste(paste("Sigma^2_",c(names(covMatList),group.names),sep="", collapse="     "), "logLik", "RSS", sep="     "))

        # store outcome as working vector
        workingY <- dat$Y
        
        # estimate variance components
        out <- .runAIREMLgaussian(	Y = workingY,
        							W = dat$W,
        							n = scan.include$n,
        							k = dat$k,
        							covMatList = covMatList,
        							m = m,
        							group.idx = group.idx,
        							g = g,
        							start = start,
        							AIREML.tol = AIREML.tol,
        							maxIter = maxIter,
        							dropZeros = dropZeros,
        							verbose = verbose)

    }else{
        # initial fit
        fit0 <- glm(dat$Y ~ -1 + dat$W, family=family)        
        eta <- fit0$linear.predictors  # W %*% beta
        mu <- family$linkinv(eta) # exp(eta)/(1 + exp(eta)) for binomial
        vmu <- family$variance(mu) # mu(1-mu) for binomial
        # inverse of g'(mu)
        gmuinv <- family$mu.eta(eta) # = vmu for canonical link
        # working vector
        workingY <- eta + (fit0$y - mu)/gmuinv
                
        Yreps <- 0
        repeat({
            Yreps <- Yreps + 1
            
            if(verbose) message(paste(paste("Sigma^2_",c(names(covMatList)),sep="", collapse="     "), "log-lik", "RSS", sep="     "))
            
            # estimate variance components
            out <- .runAIREMLother( Y = workingY, 
            						W = dat$W,
            						n = scan.include$n,
            						k = dat$k,
            						covMatList = covMatList,
            						m = m,            						
            						vmu = vmu,
            						gmuinv = gmuinv,
            						start = start,
            						AIREML.tol = AIREML.tol,
            						maxIter = maxIter,
            						dropZeros = dropZeros,           						
            						verbose = verbose)
            
            # update parameters
            if(verbose) message("Updating WorkingY Vector...")
            mu <- family$linkinv(out$eta) # exp(eta)/(1 + exp(eta)) for binomial
            vmu <- family$variance(mu) # mu(1-mu) for binomial
            # inverse of g'(mu)
            gmuinv <- family$mu.eta(out$eta) # = vmu for canonical link
            # working vector
            workingY <- out$eta + (fit0$y - mu)/gmuinv
            
            # current variance component estimate
            start <- out$varComp
            start[out$zeroFLAG] <- AIREML.tol
            
            # test for convergence
            val <- sqrt(sum((out$eta - eta)^2))
            if(verbose) message(paste("Checking for Convergence...", val, sep = "\t"))
            # update eta
            eta <- out$eta
            if(val < AIREML.tol){ 
            	break() 
            }
            if(Yreps == maxIter){
                out$converged <- FALSE
                warning("Maximum number of iterations for workingY reached without convergence!")
                break()
            }
        })
        
    }
    
    # variance component estimates 
    varComp <- out$varComp
    
    # covariance of estimates
    if(family$family == "gaussian"){
        names(varComp) <- paste("V_",c(names(covMatList),group.names),sep="")
        varCompCov <- matrix(NA, nrow=(m+g), ncol=(m+g))
        colnames(varCompCov) <- paste("V_",c(names(covMatList),group.names),sep="")
        rownames(varCompCov) <- paste("V_",c(names(covMatList),group.names),sep="")
    }else{
        names(varComp) <- paste("V_",c(names(covMatList)),sep="")
        varCompCov <- matrix(NA, nrow=m, ncol=m)
        colnames(varCompCov) <- paste("V_",c(names(covMatList)),sep="")
        rownames(varCompCov) <- paste("V_",c(names(covMatList)),sep="")
    }    
    if(dropZeros){
        varCompCov[!out$zeroFLAG, !out$zeroFLAG] <- solve(out$AI)
    }else{
        varCompCov <- solve(out$AI)
    }

    # AIC
    AIC <- 2*(dat$k + m + g) - 2*out$logLik

    # calculate conditional residuals
    residC <- as.vector(workingY - out$eta)
    
    # compute Cholesky decomposition of Covariance matrix
    cholSigmaInv <- t(chol(out$Sigma.inv))
    colnames(cholSigmaInv) <- colnames(covMatList[[1]])
    rownames(cholSigmaInv) <- rownames(covMatList[[1]])  

    # Variance Covariance of betas RSS*(Wt %*% Sigma^{-1} %*% W)^{-1}
    betaCov <- out$RSS*chol2inv(chol(crossprod(crossprod(cholSigmaInv, dat$W))))
    dimnames(betaCov) <- list(colnames(dat$W), colnames(dat$W))

    # test statistics and p-values
    SE <- sqrt(diag(betaCov))
    Stat <- (out$beta/SE)^2
    pval <- pchisq(Stat, df=1, lower.tail=FALSE)

    fixef <- data.frame(Est = out$beta, SE = SE, Stat = Stat, pval = pval)
    rownames(fixef) <- colnames(dat$W)

    family$mixedmodel <- TRUE


    return(list(varComp = varComp, 
    			varCompCov = varCompCov,
    			fixef = fixef,
    			betaCov = betaCov,
    			fitted.values = as.vector(out$fits),
    			resid.marginal = out$residM,
    			eta = as.vector(out$eta),
    			resid.conditional = residC,
    			logLikR = out$logLikR,
    			logLik = out$logLik,
    			AIC = AIC,    			
    			RSS = out$RSS,
    			workingY = as.vector(workingY),
    			model.matrix = dat$W,
    			cholSigmaInv = cholSigmaInv,
    			scanID = scan.include$value,
    			family = family,
    			converged = out$converged, 
    			zeroFLAG = out$zeroFLAG, 
    			hetResid = (g > 1) ))    
    
}



.runAIREMLgaussian <- function(Y, W, n, k, covMatList, m, group.idx, g, start, AIREML.tol, maxIter, dropZeros, verbose){
    
    # trait variance
    sigma2.p <- var(Y)
    # set convergence tolerance dependent on trait
    AIREML.tol <- AIREML.tol*sigma2.p
    val <- 2*AIREML.tol

    # initial values
    if(is.null(start)){
        sigma2.k <- rep((1/(m+1))*sigma2.p, (m+g))
    }else{
        sigma2.k <- as.vector(start)
    }
    zeroFLAG <- rep(FALSE, length(sigma2.k))
    
    reps <- 0
    repeat({
        reps <- reps+1
        
        # phenotype covariance matrix
        V <- Reduce("+", mapply("*", covMatList, sigma2.k[1:m], SIMPLIFY=FALSE))
        
        Sigma <- V
        if(g == 1){
            diag(Sigma) <- diag(Sigma) + sigma2.k[m+1]
        }else{
            diagSigma <- rep(0,n)
            for(i in 1:g){
                diagSigma[group.idx[[i]]] <- sigma2.k[m+i]
            }
            diag(Sigma) <- diag(Sigma) + diagSigma
        }
        
        # cholesky decomposition
        chol.Sigma <- chol(Sigma)
        # inverse
        Sigma.inv <- chol2inv(chol.Sigma)
        Sigma.inv_W <- crossprod(Sigma.inv,W)
        chol.Wt_Sigma.inv_W <- chol(crossprod(W, Sigma.inv_W))
        Wt_Sigma.inv_W.inv <- chol2inv(chol.Wt_Sigma.inv_W)

        # fixed effects
        beta <- crossprod(Wt_Sigma.inv_W.inv, crossprod(Sigma.inv_W,Y))
        # fitted values
        fits <- tcrossprod(W, t(beta))
        # residuals
        residM <- as.vector(Y - fits)
        # residual sum of squares
        Sigma.inv_R <- crossprod(Sigma.inv, residM)
        Rt_Sigma.inv_R <- crossprod(residM, Sigma.inv_R)        
        RSS <- as.numeric(Rt_Sigma.inv_R/(n-k))        
        # log likelihood
        logLik <- as.numeric( -0.5*n*log(2*pi*RSS) - sum(log(diag(chol.Sigma))) - 0.5*Rt_Sigma.inv_R/RSS )
        logLikR <- as.numeric( logLik + 0.5*k*log(2*pi*RSS) - sum(log(diag(chol.Wt_Sigma.inv_W))) )
        
        # print current estimates
        if(verbose) print(c(sigma2.k, logLikR, RSS))

        # check for convergence
        if(val < AIREML.tol){
        	converged <- TRUE
        	break()
        }
        if(reps > maxIter){
        	converged <- FALSE
        	warning("Maximum number of iterations reached without convergence!")
        	break()
        }
        
        # projection matrix
        P <- Sigma.inv - tcrossprod(tcrossprod(Sigma.inv_W, Wt_Sigma.inv_W.inv), Sigma.inv_W)
        PY <- crossprod(P,Y)      
        

        # run an iteration
        if(reps > 1){
            # Average Information and Scores
            AI <- matrix(NA, nrow=(m+g), ncol=(m+g))
            score <- rep(NA,(m+g))

            for(i in 1:m){
                PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
                score[i] <- -0.5*(sum(P*covMatList[[i]]) - crossprod(Y, PAPY)) # tr(PA) - YPAPY
                AI[i,i] <- 0.5*crossprod(PY, crossprod(covMatList[[i]],PAPY)) # YPAPAPY
                if((i+1) <= m){
                    for(j in (i+1):m){
                        AI[i,j] <- 0.5*crossprod(PY, crossprod(covMatList[[j]],PAPY)) # YPDPAPY
                        AI[j,i] <- AI[i,j]
                    }
                }                
                if(g == 1){
                    AI[i,(m+1)] <- 0.5*crossprod(PY, PAPY) # YPIPAPY
                    AI[(m+1),i] <- AI[i,(m+1)]
                }else{
                    for(j in 1:g){
                        AI[i,m+j] <- 0.5*crossprod(PY[group.idx[[j]]], PAPY[group.idx[[j]]]) # YP(I_group)PAPY
                        AI[m+j,i] <- AI[i,m+j]
                    }
                }
            }

            if(g == 1){
                score[m+1] <- -0.5*(sum(diag(P)) - crossprod(PY)) # tr(P) - YPIPY
                AI[(m+1),(m+1)] <- 0.5*crossprod(PY,crossprod(P,PY)) # YPIPIPY
            }else{
                for(i in 1:g){
                    PIPY <- crossprod(P[group.idx[[i]], ],PY[group.idx[[i]]])
                    score[m+i] <- -0.5*(sum(diag(P)[group.idx[[i]]]) - crossprod(PY[group.idx[[i]]])) # tr(P(I_group)) - YP(I_group)PY
                    AI[m+i,m+i] <- 0.5*crossprod(PY[group.idx[[i]]], PIPY[group.idx[[i]]]) # YP(I_group)P(I_group)PY
                    if((i+1) <= g){
                        for(j in (i+1):g){
                            AI[m+i,m+j] <- 0.5*crossprod(PY[group.idx[[j]]], PIPY[group.idx[[j]]]) # YP(I_group2)P(I_group)PY
                            AI[m+j,m+i] <- AI[m+i,m+j]
                        }
                    }
                }
            }
            
            if(dropZeros){
                # remove Zero terms
                AI <- AI[!zeroFLAG,!zeroFLAG]
                score <- score[!zeroFLAG]
            }
            
            # update
            AIinvScore <- solve(AI, score)
            
            if(dropZeros){
                sigma2.kplus1[!zeroFLAG] <- sigma2.k[!zeroFLAG] + AIinvScore
                sigma2.kplus1[zeroFLAG] <- 0
            }else{
                sigma2.kplus1 <- sigma2.k + AIinvScore
                # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
                sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 
            }
            
            # step-halving if step too far
            tau <- 1
            while(!all(sigma2.kplus1 >= 0)){
                tau <- 0.5*tau
                if(dropZeros){
                    sigma2.kplus1[!zeroFLAG] <- sigma2.k[!zeroFLAG] + tau*AIinvScore
                    sigma2.kplus1[zeroFLAG] <- 0
                }else{
                    sigma2.kplus1 <- sigma2.k + tau*AIinvScore
                    sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 
                }
            }

            # check for zeros
            zeroFLAG <- sigma2.kplus1 < AIREML.tol # which elements have converged to "0"
            sigma2.kplus1[zeroFLAG] <- 0 # set these to 0
            
            # test for convergence
            val <- sqrt(sum((sigma2.kplus1 - sigma2.k)^2))
            # update estimates
            sigma2.k <- sigma2.kplus1
            
            
        }else{
            # EM step
            sigma2.kplus1 <- rep(NA,(m+g))
            for(i in 1:m){
                PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
                sigma2.kplus1[i] <- (1/n)*(sigma2.k[i]^2*crossprod(Y,PAPY) + n*sigma2.k[i] - sigma2.k[i]^2*sum(P*covMatList[[i]]))
            }
            if(g == 1){
                sigma2.kplus1[m+1] <- (1/n)*(sigma2.k[m+1]^2*crossprod(PY) + n*sigma2.k[m+1] - sigma2.k[m+1]^2*sum(diag(P)))
            }else{
                for(i in 1:g){
                    sigma2.kplus1[m+i] <- (1/n)*(sigma2.k[m+i]^2*crossprod(PY[group.idx[[i]]]) + n*sigma2.k[m+i] - sigma2.k[m+i]^2*sum(diag(P)[group.idx[[i]]]))
                }
            }
            # update estimates
            sigma2.k <- sigma2.kplus1
        }
        
    })
    
    # linear predictor
    eta <- fits + crossprod(V, Sigma.inv_R) # X\beta + Zb
    
    return(list(varComp = sigma2.k, 
    			AI = AI, 
    			converged = converged, 
    			zeroFLAG = zeroFLAG,
    			Sigma.inv = Sigma.inv,
    			beta = beta,
    			residM = residM,
    			fits = fits,
    			eta = eta,
    			logLikR = logLikR,
    			logLik = logLik,
    			RSS=RSS))    
}



.runAIREMLother <- function(Y, W, n, k, covMatList, m, vmu, gmuinv, start, AIREML.tol, maxIter, dropZeros, verbose){
    
    val <- 2*AIREML.tol

    # initial values for variance components
    if(is.null(start)){
        sigma2.k <- rep(10*sqrt(AIREML.tol), m)
    }else{
        sigma2.k <- as.vector(start)
    }
    zeroFLAG <- rep(FALSE, length(sigma2.k))
    
    reps <- 0
    repeat({
        reps <- reps+1
                
        # phenotype covariance matrix
        V <- Reduce("+", mapply("*", covMatList, sigma2.k[1:m], SIMPLIFY=FALSE))
        Sigma <- V + diag(as.vector(vmu)/as.vector(gmuinv)^2)

        # cholesky decomposition
        chol.Sigma <- chol(Sigma)
        # inverse
        Sigma.inv <- chol2inv(chol.Sigma)
        Sigma.inv_W <- crossprod(Sigma.inv,W)
        chol.Wt_Sigma.inv_W <- chol(crossprod(W, Sigma.inv_W))
        Wt_Sigma.inv_W.inv <- chol2inv(chol.Wt_Sigma.inv_W)

        # fixed effects
        beta <- crossprod(Wt_Sigma.inv_W.inv, crossprod(Sigma.inv_W,Y))
        # fitted values
        fits <- tcrossprod(W, t(beta))
        # residuals
        residM <- as.vector(Y - fits)
        # residual sum of squares
        Sigma.inv_R <- crossprod(Sigma.inv, residM)
        Rt_Sigma.inv_R <- crossprod(residM, Sigma.inv_R)        
        RSS <- as.numeric(Rt_Sigma.inv_R/(n-k))        
        # log likelihood
        logLik <- as.numeric( -0.5*n*log(2*pi*RSS) - sum(log(diag(chol.Sigma))) - 0.5*Rt_Sigma.inv_R/RSS )
        logLikR <- as.numeric( logLik + 0.5*k*log(2*pi*RSS) - sum(log(diag(chol.Wt_Sigma.inv_W))) )
        
        # print current estimates
        if(verbose) print(c(sigma2.k, logLikR, RSS))

        # check for convergence
        if(val < AIREML.tol){
        	converged <- TRUE
        	break()
        }
        if(reps > maxIter){
        	converged <- FALSE
        	warning("Maximum number of iterations reached without convergence!")
        	break()
        }

        # projection matrix
        P <- Sigma.inv - tcrossprod(tcrossprod(Sigma.inv_W, Wt_Sigma.inv_W.inv), Sigma.inv_W)
        PY <- crossprod(P,Y)   


        # run an iteration
        if(reps > 1){
            # Average Information and Scores
            AI <- matrix(NA, nrow=m, ncol=m)
            score <- rep(NA,m)

            for(i in 1:m){
                PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
                score[i] <- -0.5*(sum(P*covMatList[[i]]) - crossprod(Y, PAPY))
                AI[i,i] <- 0.5*crossprod(PY, crossprod(covMatList[[i]],PAPY)) # YPAPAPY
                if((i+1) <= m){
                    for(j in (i+1):m){
                        AI[i,j] <- 0.5*crossprod(PY, crossprod(covMatList[[j]],PAPY)) # YPDPAPY
                        AI[j,i] <- AI[i,j]
                    }
                }
            }
            
            if(dropZeros){
                # remove Zero terms
                AI <- AI[!zeroFLAG,!zeroFLAG]
                score <- score[!zeroFLAG]
            }
            
            # update
            AIinvScore <- solve(AI, score)
            
            if(dropZeros){
                sigma2.kplus1[!zeroFLAG] <- sigma2.k[!zeroFLAG] + AIinvScore
                sigma2.kplus1[zeroFLAG] <- 0
            }else{
                sigma2.kplus1 <- sigma2.k + AIinvScore
                # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
                sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 
            }
            
            # step-halving if step too far
            tau <- 1
            while(!all(sigma2.kplus1 >= 0)){
                tau <- 0.5*tau
                if(dropZeros){
                    sigma2.kplus1[!zeroFLAG] <- sigma2.k[!zeroFLAG] + tau*AIinvScore
                    sigma2.kplus1[zeroFLAG] <- 0
                }else{
                    sigma2.kplus1 <- sigma2.k + tau*AIinvScore
                    sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 
                }
            }

            # check for zeros
            zeroFLAG <- sigma2.kplus1 < AIREML.tol # which elements have converged to "0"
            sigma2.kplus1[zeroFLAG] <- 0 # set these to 0
            
            # test for convergence
            val <- sqrt(sum((sigma2.kplus1 - sigma2.k)^2))
            # update estimates
            sigma2.k <- sigma2.kplus1
            
        }else{
            # EM step
            sigma2.kplus1 <- rep(NA,m)
            for(i in 1:m){
                PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
                sigma2.kplus1[i] <- (1/n)*((sigma2.k[i])^2*crossprod(Y,PAPY) + (n*sigma2.k[i] - (sigma2.k[i])^2*sum(P*covMatList[[i]])))
            }
            # update estimates
            sigma2.k <- sigma2.kplus1
        }
        
    })
    
    # linear predictor
    eta <- fits + crossprod(V, Sigma.inv_R) # X\beta + Zb

    return(list(varComp = sigma2.k, 
    			AI = AI, 
    			converged = converged, 
    			zeroFLAG = zeroFLAG,
    			Sigma.inv = Sigma.inv,
    			beta = beta,
    			residM = residM,
    			fits = fits,
    			eta = eta,
    			logLikR = logLikR,
    			logLik = logLik,
    			RSS=RSS)) 
}


