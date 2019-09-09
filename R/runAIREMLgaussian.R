.runAIREMLgaussian <- function(Y, X, start, covMatList, group.idx, AIREML.tol, 
                                drop.zeros, max.iter, EM.iter, verbose){

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
    if(verbose) message(paste(paste("Sigma^2_",c(names(covMatList)),sep="", collapse="     "), "log-lik", "RSS", sep="     "))

    reps <- 0
    repeat({
        reps <- reps+1
        
        ### compute sigma quantities
        sq <- .computeSigmaQuantities(varComp = sigma2.k, covMatList = covMatList, group.idx = group.idx)
        ### compute likelihood quantities
        lq <- .calcLikelihoodQuantities(Y = Y, X = X, Sigma.inv = sq$Sigma.inv, cholSigma.diag = sq$cholSigma.diag)

        # print current estimates
        if(verbose) print(c(sigma2.k, lq$logLikR, lq$RSS))
                
        if(reps > EM.iter){
            # Average Information and Scores
            AI <- matrix(NA, nrow=(m+g), ncol=(m+g))
            score <- rep(NA,(m+g))
            covMats.score.AI <- .calcAIcovMats(PY = lq$PY, covMatList = covMatList,
                                               Sigma.inv = sq$Sigma.inv, Sigma.inv_X = lq$Sigma.inv_X, Xt_Sigma.inv_X.inv = lq$Xt_Sigma.inv_X.inv)
            AI[1:m, 1:m] <- covMats.score.AI$AI
            score[1:m]  <- covMats.score.AI$score
            het.vars.score.AI <- .calcAIhetvars(PY = lq$PY, group.idx = group.idx,
                                                Sigma.inv = sq$Sigma.inv, Sigma.inv_X = lq$Sigma.inv_X, Xt_Sigma.inv_X.inv = lq$Xt_Sigma.inv_X.inv)
            score[(m + 1):(m + g)]  <- het.vars.score.AI$score
            AI[(m + 1):(m + g),(m+1):(m + g)]  <- het.vars.score.AI$AI
            
            ### take care of "off diagonal" (terms for covariance between variance components corresponding to 
            ### the random effects and the residuals variances) 
            AI.off <- .calcAIcovMatsResids(PY = lq$PY, covMatList = covMatList, group.idx = group.idx,
                                           Sigma.inv = sq$Sigma.inv, Sigma.inv_X = lq$Sigma.inv_X, Xt_Sigma.inv_X.inv = lq$Xt_Sigma.inv_X.inv)
            AI[1:m, (m + 1):(m + g)] <- AI.off
            AI[(m + 1):(m + g), 1:m] <- t(AI.off)
            
            if(drop.zeros){
                # remove Zero terms
                AI <- AI[!zeroFLAG,!zeroFLAG]
                score <- score[!zeroFLAG]
            }
            
            # update
            AIinvScore <- solve(AI, score)
            
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
                # PAPY <- sq$Sigma.inv %*% crossprod(covMatList[[i]],lq$PY) - tcrossprod(tcrossprod(lq$Sigma.inv_X, lq$Xt_Sigma.inv_X.inv), t(crossprod(covMatList[[i]],lq$PY)) %*% lq$Sigma.inv_X)
                trPA.part1 <- sum( sq$Sigma.inv * covMatList[[i]] )
                trPA.part2 <- sum(diag( (crossprod( lq$Sigma.inv_X, covMatList[[i]]) %*% lq$Sigma.inv_X) %*% lq$Xt_Sigma.inv_X.inv ))
                trPA <-  trPA.part1 - trPA.part2
                
                APY <- crossprod(covMatList[[i]],lq$PY)
                sigma2.kplus1[i] <- as.numeric((1/n)*(sigma2.k[i]^2*crossprod(lq$PY,APY) + n*sigma2.k[i] - sigma2.k[i]^2*trPA ))
                # sigma2.kplus1[i] <- as.numeric((1/n)*(sigma2.k[i]^2*crossprod(Y,PAPY) + n*sigma2.k[i] - sigma2.k[i]^2*trPA )) 
            }
            if(g == 1){    
                trP.part1 <- sum(diag( sq$Sigma.inv ))
                trP.part2 <- sum(diag( crossprod( lq$Sigma.inv_X) %*% lq$Xt_Sigma.inv_X.inv ))
                trP <-  trP.part1 - trP.part2
                
                sigma2.kplus1[m+1] <- as.numeric((1/n)*(sigma2.k[m+1]^2*crossprod(lq$PY) + n*sigma2.k[m+1] - sigma2.k[m+1]^2*trP ))
            }else{
                for(i in 1:g){
                    covMati <- as.numeric( 1:n %in% group.idx[[i]] )
                    trPi.part1 <- sum(diag(sq$Sigma.inv)[ group.idx[[i]] ] )
                    trPi.part2 <- sum(diag( crossprod(crossprod(lq$Sigma.inv_X*covMati, lq$Sigma.inv_X), lq$Xt_Sigma.inv_X.inv )))
                    trPi <- trPi.part1 - trPi.part2
                    sigma2.kplus1[m+i] <- as.numeric((1/n)*(sigma2.k[m+i]^2*crossprod(lq$PY[group.idx[[i]]]) + n*sigma2.k[m+i] - sigma2.k[m+i]^2*trPi ))
                }
            }
        }

        ### check for convergence
        # val <- sqrt(sum((sigma2.kplus1 - sigma2.k)^2))
        if((reps > EM.iter) & (max(abs(sigma2.kplus1 - sigma2.k)) < AIREML.tol)){
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

    # linear predictor
    eta <- as.numeric(lq$fits + crossprod(sq$Vre, lq$PY)) # X\beta + Zb
    
    return(list(varComp = sigma2.k, AI = AI, converged = converged, zeroFLAG = zeroFLAG, niter = reps,
                Sigma.inv = sq$Sigma.inv, beta = lq$beta, residM = lq$residM, fits = lq$fits, eta = eta, 
                logLikR = lq$logLikR, logLik = lq$logLik, RSS = lq$RSS))  
}





