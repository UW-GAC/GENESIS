.runAIREMLother <- function(Y, X, start, covMatList, AIREML.tol, drop.zeros, max.iter, EM.iter, verbose, vmu, gmuinv){
    
    m <- length(covMatList)
    n <- length(Y)
    
    # initial values for variance components
    if(is.null(start)){
        sigma2.k <- rep(10*sqrt(AIREML.tol), m)
    }else{
        sigma2.k <- as.vector(start)
    }
    
    reps <- 0
    repeat({
        reps <- reps+1
        
        zeroFLAG <- (sigma2.k < AIREML.tol) # which elements have converged to "0"
        sigma2.k[zeroFLAG] <- 0 # set these to 0
        
        if (sum(zeroFLAG) == m)  return(list(allZero = TRUE))
        
        sq <- .computeSigmaQuantities(varComp = sigma2.k, covMatList = covMatList, vmu = vmu, gmuinv = gmuinv )     
        lq <- .calcLikelihoodQuantities(Y, X, sq$Sigma.inv, diag(sq$cholSigma))


        
        # print current estimates
        if(verbose) print(c(sigma2.k, lq$logLikR, lq$RSS))

        if(reps > EM.iter){
            # Average Information and Scores
            ### more arguments
            covMats.score.AI <- .calcAIcovMats(Y, lq$PY, covMatList,
                                               Sigma.inv = sq$Sigma.inv, Sigma.inv_X = lq$Sigma.inv_X, Xt_Sigma.inv_X.inv = lq$Xt_Sigma.inv_X.inv)
            AI <- covMats.score.AI$AI
            score <- covMats.score.AI$score
            
            if(drop.zeros){  ## here need to exit if all terms were zero!!
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
                sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
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
                    sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
                }
            }
            
            # test for convergence
            stat <- sqrt(sum((sigma2.kplus1 - sigma2.k)^2))
            sigma2.k <- sigma2.kplus1
            if(stat < AIREML.tol){
                converged <- TRUE
                break()
            }
            if(reps == max.iter){
                converged <- FALSE
                warning("Maximum number of iterations reached without convergence!")
                break()
            }
            
        }else{
            # EM step
            sigma2.kplus1 <- rep(NA,m)
            for(i in 1:m){
                ### PAPY <- crossprod(lq$P,crossprod(covMatList[[i]],lq$PY))
                ### sigma2.kplus1[i] <- (1/n)*((sigma2.k[i])^2*crossprod(Y,lq$PAPY) + (n*sigma2.k[i] - (sigma2.k[i])^2*sum(lq$P*covMatList[[i]])))
                PAPY <- sq$Sigma.inv %*% crossprod(covMatList[[i]],lq$PY) - tcrossprod(tcrossprod(lq$Sigma.inv_X, lq$Xt_Sigma.inv_X.inv), t(crossprod(covMatList[[i]],lq$PY)) %*% lq$Sigma.inv_X)	  
                trPA.part1 <- sum( sq$Sigma.inv * covMatList[[i]] )
                trPA.part2 <- sum(diag( (crossprod( lq$Sigma.inv_X, covMatList[[i]]) %*% lq$Sigma.inv_X) %*% lq$Xt_Sigma.inv_X.inv ))
                trPA <-  trPA.part1 - trPA.part2
                
                sigma2.kplus1[i] <- as.numeric((1/n)*(sigma2.k[i]^2*crossprod(Y,PAPY) + n*sigma2.k[i] - sigma2.k[i]^2*trPA ))
                
            }
            sigma2.k <- sigma2.kplus1
        }
        
    })
    
    
    # linear predictor
    VinvR <- crossprod(sq$Sigma.inv, lq$residM)
    eta <- as.numeric(lq$fits + crossprod(sq$Vre, VinvR)) # X\beta + Zb
    
    return(list(allZero = FALSE, varComp = sigma2.k, AI = AI, converged = converged, zeroFLAG = zeroFLAG, beta = lq$beta, residM = lq$residM, 
                eta = eta, logLikR = lq$logLikR, logLik = lq$logLik, RSS = lq$RSS, fits = lq$fits, iter = reps))

}


