.runAIREMLgaussianPCG <- function(Y, X, start, covMatList, group.idx, AIREML.tol, 
                                    drop.zeros, max.iter, EM.iter, verbose, Method = "EST"){
    
    # initial values
    m <- length(covMatList)
    g <- length(group.idx)
    n <- length(Y)
    sigma2.p <- drop(var(Y))
    AIREML.tol <- AIREML.tol*sigma2.p  # set convergence tolerance dependent on trait
    if(is.null(start)){
        sigma2.k <- rep((1/(m+1))*sigma2.p, (m+g))
    }else{
        sigma2.k <- as.vector(start)
    }
    sigma2.kplus1 <- rep(NA, length(sigma2.k))
    zeroFLAG <- rep(FALSE, length(sigma2.k))

    if(verbose) message("Computing Variance Component Estimates...")
    if(verbose) message(paste("Sigma^2_",c(names(covMatList)),sep="", collapse="     "))

    reps <- 0
    repeat({
        reps <- reps+1
        
        ### compute sigma quantities
        sq <- .computeSigma(varComp = sigma2.k, covMatList = covMatList, group.idx = group.idx)
        Sigma.inv_X <- .pcgM(sq$Sigma, X)
        Sigma <- sq$Sigma
        PY <- .LeftMP(Y, X, Sigma, Sigma.inv_X)
        
        # lq <- .calcLikelihoodQuantitiesPCG(Y, X, Sigma, Sigma.inv_X, PPY=(g==1))
        
        # print current estimates
        if(verbose) print(sigma2.k)
        
        if(reps > EM.iter){
            # Average Information and Scores
            AI <- matrix(NA, nrow=(m+g), ncol=(m+g))
            score <- rep(NA,(m+g))
            covMats.score.AI <- .calcAIcovMatsPCG(X, PY, covMatList, Sigma, Sigma.inv_X, Method = Method)
            AI[1:m, 1:m] <- covMats.score.AI$AI
            score[1:m]  <- covMats.score.AI$score
            het.vars.score.AI <- .calcAIhetvarsPCG(X, PY, covMatList, Sigma, Sigma.inv_X,group.idx, Method = Method)
            score[(m + 1):(m + g)]  <- het.vars.score.AI$score
            AI[(m + 1):(m + g),(m+1):(m + g)]  <- het.vars.score.AI$AI
            
            ### take care of "off diagonal" (terms for covariance between variance components corresponding to 
            ### the covariance matriecs, and the residuals variances) 
            AI.off <- .calcAIcovMatsResidsPCG(X, PY, covMatList, Sigma, Sigma.inv_X,group.idx)
            AI[1:m, (m + 1):(m + g)] <- AI.off
            AI[(m + 1):(m + g),1:m ] <- t(AI.off)
            
            if(drop.zeros){
                # remove Zero terms
                AI <- AI[!zeroFLAG,!zeroFLAG]
                score <- score[!zeroFLAG]
            }
            
            # update
            AIinvScore <- solve(AI, score)
            # print("score is")
            # print(score)
            # print("AI INV is")
            # print(solve(AI))
            if(drop.zeros){
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
                if(drop.zeros){
                    sigma2.kplus1[!zeroFLAG] <- sigma2.k[!zeroFLAG] + tau*AIinvScore
                    sigma2.kplus1[zeroFLAG] <- 0
                }else{
                    sigma2.kplus1 <- sigma2.k + tau*AIinvScore
                    # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
                    sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 
                }
            }
            
        }else{
            # EM step
            for(i in 1:m){
                ### PAPY <- crossprod(lq$P,crossprod(covMatList[[i]],lq$PY))
                APY <- crossprod(covMatList[[i]],PY)
                # PAPY <- .LeftMP(APY,X, Sigma, Sigma.inv_X)
                trPA <-  .LeftMPtrEst(covMatList[[i]], X, Sigma, Sigma.inv_X,Method = Method)
                ### sigma2.kplus1[i] <- (1/n)*(sigma2.k[i]^2*crossprod(Y,PAPY) + n*sigma2.k[i] - sigma2.k[i]^2*sum(lq$P*covMatList[[i]]))
                sigma2.kplus1[i] <- as.numeric((1/n)*(sigma2.k[i]^2*crossprod(PY,APY) + n*sigma2.k[i] - sigma2.k[i]^2*trPA ))
                
            }
            if(g == 1){
                ### sigma2.kplus1[m+1] <- (1/n)*(sigma2.k[m+1]^2*crossprod(lq$PY) + n*sigma2.k[m+1] - sigma2.k[m+1]^2*sum(diag(lq$P)))
                
                I <- diag(n)
                trP <-  .LeftMPtrEst(I, X, Sigma, Sigma.inv_X,Method = Method)
                
                sigma2.kplus1[m+1] <- as.numeric((1/n)*(sigma2.k[m+1]^2*crossprod(PY) + n*sigma2.k[m+1] - sigma2.k[m+1]^2*trP ))
            }else{
                for(i in 1:g){
                    ### sigma2.kplus1[m+i] <- (1/n)*(sigma2.k[m+i]^2*crossprod(lq$PY[group.idx[[i]]]) + n*sigma2.k[m+i] - sigma2.k[m+i]^2*sum(diag(lq$P)[group.idx[[i]]]))
                    ## covMati <- Diagonal( x=as.numeric( 1:n %in% group.idx[[i]] ) )
                    ## trPi.part1 <- sum(diag(sq$Sigma.inv)[ group.idx[[i]] ] )
                    ## trPi.part2 <- sum(diag( (crossprod( lq$Sigma.inv_X, covMati) %*% lq$Sigma.inv_X) %*% lq$Xt_Sigma.inv_X.inv ))
                    I <- as.numeric( 1:n %in% group.idx[[i]] )
                    trPi <- .LeftMPtrEst(I, X, Sigma, Sigma.inv_X,Method = Method)
                    sigma2.kplus1[m+i] <- as.numeric((1/n)*(sigma2.k[m+i]^2*crossprod(PY[group.idx[[i]]]) + n*sigma2.k[m+i] - sigma2.k[m+i]^2*trPi ))
                }
            }
        }

        ### check for convergence
        if((reps > EM.iter) & (max(abs(sigma2.k - sigma2.kplus1)/(abs(sigma2.k) + abs(sigma2.kplus1) + 0.02)) < 0.02)){
            converged <- TRUE
            (break)()
        }else{
            # check if exceeded the number of iterations
            if(reps == max.iter){
                converged <- FALSE
                warning("Maximum number of iterations reached without convergence!")
                (break)()
            }else{
                # check which parameters have converged to "0"
                zeroFLAG <- sigma2.kplus1 < AIREML.tol
                sigma2.kplus1[zeroFLAG] <- 0
                # update estimates
                sigma2.k <- sigma2.kplus1
            }
        }
    })

    ### some of the things normally in calcLQ ###

    # get beta estimates
    Xt_Sigma.inv_X <- crossprod(X, Sigma.inv_X)
    # fix issue with not recognizing the matrix as symmetric
    Xt_Sigma.inv_X <- (Xt_Sigma.inv_X + t(Xt_Sigma.inv_X))/2
    chol.Xt_Sigma.inv_X <- chol(Xt_Sigma.inv_X)
    Xt_Sigma.inv_X.inv <- chol2inv(chol.Xt_Sigma.inv_X)
    beta <- crossprod(Xt_Sigma.inv_X.inv, crossprod(Sigma.inv_X, Y))

    ## calculate the mean of the outcomes
    fits <- tcrossprod(X, t(beta))
    
    # linear predictor
    eta <- as.numeric(fits + crossprod(sq$Vre, PY)) # X\beta + Zb
    
    # obtain marginal residuals
    residM <- as.vector(Y - fits)

    # compute RSS
    YPY <- crossprod(Y, PY)
    RSS <- as.numeric(YPY/(n-ncol(X)))

    # get Sigma.inv
    Sigma.inv <- chol2inv(chol(sq$Sigma))
    
    # lq <- .calcLikelihoodQuantitiesPCG(Y, X, Sigma, Sigma.inv_X, PPY=(g==1))
    
    return(list(varComp = sigma2.k, AI = AI, converged = converged, zeroFLAG = zeroFLAG, niter = reps, 
                Sigma.inv = Sigma.inv, beta = beta, residM = residM, fits = fits, eta = eta, 
                logLikR = NULL, logLik = NULL, RSS = RSS))

                ### do we need Sigma.inv here? ### for .nullModOutMM?    
}




