.runAIREMLgaussianPCGavg <- function(Y, X, start, covMatList, group.idx, AIREML.tol, drop.zeros, max.iter, verbose, Method = "EST"){
    # initial values
    m <- length(covMatList)
    #print(m)
    g <- length(group.idx)
    #print(g)
    n <- length(Y)
    sigma2.p <- drop(var(Y))
    AIREML.tol <- AIREML.tol*sigma2.p  # set convergence tolerance dependent on trait
    val <- 2 * AIREML.tol
    if(is.null(start)){
        sigma2.k <- rep((1/(m+1))*sigma2.p, (m+g))
    }else{
        sigma2.k <- as.vector(start)
    }

    reps <- 0
    
    repeat({
        reps <- reps+1
        
        zeroFLAG <- sigma2.k < AIREML.tol # which elements have converged to "0"
        sigma2.k[zeroFLAG] <- 0 # set these to 0
        
        Sigma <- .computeSigma(varComp = sigma2.k, covMatList = covMatList, group.idx = group.idx)
        Sigma.inv_X <- .pcgM(Sigma, X)
        PY <- .LeftMP(Y, X, Sigma, Sigma.inv_X)
        
        # print current estimates
        if(verbose) print(c(sigma2.k))
        
        REMLPCGavg[reps] <<- .calREML(Y, X, Sigma, Sigma.inv_X, PY)
        ## check for convergence
        if (val < AIREML.tol) {
            converged <- TRUE
            (break)()
        }
        
        ## check if exceeded the number of iterations
        if (reps > max.iter) {
            converged <- FALSE
            warning("Maximum number of iterations reached without convergence!")
            (break)()
        }
        
        
        if(reps > 1){
            # Average Information and Scores
            AI <- matrix(NA, nrow=(m+g), ncol=(m+g))
            score <- rep(NA,(m+g))
            covMats.score.AI <- .calcAIcovMatsPCG(Y, X, PY, covMatList, Sigma, Sigma.inv_X, Method = Method)
            AI[1:m, 1:m] <- covMats.score.AI$AI
            score[1:m]  <- covMats.score.AI$score
            het.vars.score.AI <- .calcAIhetvarsPCG(Y, X, PY, covMatList, Sigma, Sigma.inv_X,group.idx, Method = Method)
            score[(m + 1):(m + g)]  <- het.vars.score.AI$score
            AI[(m + 1):(m + g),(m+1):(m + g)]  <- het.vars.score.AI$AI
            
            ### take care of "off diagonal" (terms for covariance between variance components corresponding to 
            ### the covariance matriecs, and the residuals variances) 
            AI.off <- .calcAIcovMatsResidsPCG(Y, X, PY, covMatList, Sigma, Sigma.inv_X,group.idx)
            AI[1:m, (m + 1):(m + g)] <- AI.off
            AI[(m + 1):(m + g),1:m ] <- t(AI.off)
            
            if(drop.zeros){
                # remove Zero terms
                AI <- AI[!zeroFLAG,!zeroFLAG]
                score <- score[!zeroFLAG]
            }
            
            # update
            AIinvScore <- solve(AI, score)
            print("score is")
            print(score)
            print("AI INV is")
            print(solve(AI))
            if(drop.zeros){
                sigma2.kplus1[!zeroFLAG] <- sigma2.k[!zeroFLAG] + AIinvScore/reps
                sigma2.kplus1[zeroFLAG] <- 0
            }else{
                sigma2.kplus1 <- sigma2.k + AIinvScore/reps
                sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
            }
            
            # step-halving if step too far
            # tau <- 1
            # while(!all(sigma2.kplus1 >= 0)){
            #     tau <- 0.5*tau
            #     if(drop.zeros){
            #         sigma2.kplus1[!zeroFLAG] <- sigma2.k[!zeroFLAG]+ tau*AIinvScore/reps
            #         sigma2.kplus1[zeroFLAG] <- 0
            #     }else{
            #         sigma2.kplus1 <- sigma2.k + tau*AIinvScore/reps
            #         sigma2.kplus1[zeroFLAG & sigma2.kplus1 < AIREML.tol] <- 0 # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
            #     }
            # }
            
            zeroFLAG <- sigma2.kplus1 < AIREML.tol
            sigma2.kplus1[zeroFLAG] <- 0
            val <- sqrt(sum((sigma2.kplus1 - sigma2.k)^2))

            # update estimates
            sigma2.k <- sigma2.kplus1
            print(sigma2.k)
        }else{
            # EM step
            sigma2.kplus1 <- rep(NA,(m+g))
            for(i in 1:m){
                ### PAPY <- crossprod(lq$P,crossprod(covMatList[[i]],lq$PY))
                APY <- crossprod(covMatList[[i]],PY)
                PAPY <- .LeftMP(APY,X, Sigma, Sigma.inv_X)
                trPA <-  .LeftMPtrEst(covMatList[[i]], X, Sigma, Sigma.inv_X,Method = Method)
                
                ### sigma2.kplus1[i] <- (1/n)*(sigma2.k[i]^2*crossprod(Y,PAPY) + n*sigma2.k[i] - sigma2.k[i]^2*sum(lq$P*covMatList[[i]]))
                sigma2.kplus1[i] <- as.numeric((1/n)*(sigma2.k[i]^2*crossprod(Y,PAPY) + n*sigma2.k[i] - sigma2.k[i]^2*trPA ))
                
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
                    sigma2.kplus1[m+i] <- as.numeric((1/n)*(sigma2.k[m+i]^2*crossprod(lq$PY[group.idx[[i]]]) + n*sigma2.k[m+i] - sigma2.k[m+i]^2*trPi ))
                }
            }
            sigma2.k <- sigma2.kplus1
        }
        
    })
    
    # linear predictor
    #eta <- as.numeric(lq$fits + crossprod(sq$Vre, lq$Sigma.inv_R)) # X\beta + Zb
    print("reach the end of runAIREMLgaussianPCGavg")
    return(list(varComp = sigma2.k, AI = AI, converged = converged, zeroFLAG = zeroFLAG))
    
}





