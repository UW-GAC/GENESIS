

.calcAIcovMatsMulti <- function(PY, covMatList, Sigma.inv, Sigma.inv_X, Xt_Sigma.inv_X.inv, nsamp, np, kinMatSmall){
    
    # There are np partial derivatives for the AI matrix calculation. They are all block diagonal, with mostly 0’s and either identity or kinship matrix in given quadrants/blocks.
    # First set of partial derivatives are the trait-specific resid terms. 
    # These are block-diagonal matrices with the identity on the diagonal block of each given trait, 
    # ie. trait 1 has identity in block 1,1, trait 2 had identity in block 2,2, …
    
    m <- length(covMatList)
    AI <- matrix(NA, nrow =  np, ncol = np)
    score <- rep(NA, np)
    
    zeroMat <- Matrix(0,nrow=nsamp,ncol=nsamp)
    identMat <- zeroMat+diag(nsamp)
    I1 <- bdiag(identMat,zeroMat)
    I2 <- bdiag(zeroMat,identMat)
    
    # Second set of partial derivatives are the trait-specific genotype covar terms, 
    # ie. the kinship matrix. 
    # These are block diagonal matrices with the kinship in the block of each given trait, 
    # ie. trait 1 has kinship in block 1,1, trait 2 has kinship in block 2,2, …
    
    G1 <- bdiag(kinMatSmall,zeroMat)
    G2 <- bdiag(zeroMat,kinMatSmall)
    
    # Finally, the cross-trait genotype covar term. This has kinship on the off-diagonals and 0’s elsewhere.
    tmpMat <- matrix(c(0,1,1,0),nrow=2)
    C <- kronecker(tmpMat,kinMatSmall)
    
    # Now we have all the AI elements calculated. Next we need to multiply these by the projection matrix P
    I1PY <- crossprod(I1,PY)
    I2PY <- crossprod(I2,PY)
    G1PY <- crossprod(G1,PY)
    G2PY <- crossprod(G2,PY)
    CPY <- crossprod(C,PY)
    
    PI1PY <- crossprod(Sigma.inv, I1PY) - crossprod(tcrossprod(Xt_Sigma.inv_X.inv, Sigma.inv_X), crossprod(Sigma.inv_X, I1PY))
    PI2PY <- crossprod(Sigma.inv, I2PY) - crossprod(tcrossprod(Xt_Sigma.inv_X.inv, Sigma.inv_X), crossprod(Sigma.inv_X, I2PY))
    PG1PY <- crossprod(Sigma.inv, G1PY) - crossprod(tcrossprod(Xt_Sigma.inv_X.inv, Sigma.inv_X), crossprod(Sigma.inv_X, G1PY))
    PG2PY <- crossprod(Sigma.inv, G2PY) - crossprod(tcrossprod(Xt_Sigma.inv_X.inv, Sigma.inv_X), crossprod(Sigma.inv_X, G2PY))
    PCPY <- crossprod(Sigma.inv, CPY) - crossprod(tcrossprod(Xt_Sigma.inv_X.inv, Sigma.inv_X), crossprod(Sigma.inv_X, CPY))
    
    trPI1.part1 <- sum( Sigma.inv * I1 )
    trPI1.part2 <- sum(diag( (crossprod( Sigma.inv_X, I1) %*% Sigma.inv_X) %*% Xt_Sigma.inv_X.inv ))
    trPI1 <-  trPI1.part1 - trPI1.part2
    
    trPI2.part1 <- sum( Sigma.inv * I2 )
    trPI2.part2 <- sum(diag( (crossprod( Sigma.inv_X, I2) %*% Sigma.inv_X) %*% Xt_Sigma.inv_X.inv ))
    trPI2 <-  trPI2.part1 - trPI2.part2
    
    trPG1.part1 <- sum( Sigma.inv * G1 )
    trPG1.part2 <- sum(diag( (crossprod( Sigma.inv_X, G1) %*% Sigma.inv_X) %*% Xt_Sigma.inv_X.inv ))
    trPG1 <-  trPG1.part1 - trPG1.part2
    
    trPG2.part1 <- sum( Sigma.inv * G2 )
    trPG2.part2 <- sum(diag( (crossprod( Sigma.inv_X, G2) %*% Sigma.inv_X) %*% Xt_Sigma.inv_X.inv ))
    trPG2 <-  trPG2.part1 - trPG2.part2
    
    trPC.part1 <- sum( Sigma.inv * C )
    trPC.part2 <- sum(diag( (crossprod( Sigma.inv_X, C) %*% Sigma.inv_X) %*% Xt_Sigma.inv_X.inv ))
    trPC <-  trPC.part1 - trPC.part2
    
    allTR <- c(trPI1,trPI2,trPG1,trPG2,trPC)
    
    APY <- list(I1PY, I2PY, G1PY, G2PY, CPY)
    PAPY <- list(PI1PY, PI2PY, PG1PY, PG2PY, PCPY)
    
    # these are the commands for the single trait analysis calculation of score and AI
    # score[i] <- as.numeric(-0.5*(trPA - crossprod(PY, APY))) # tr(PA) - YPAPY
    # AI[i,i] <- as.numeric(0.5*crossprod(APY, PAPY)) # YPAPAPY
    
    for(i in 1:length(score)){
        score[i] <- (-0.5*(allTR[i] - crossprod(PY, APY[[i]])))
    }
    for(i in 1:length(score)){ # this is the row index
        for(j in 1:i){ # this is the column index
            AI[i,j] <- as.numeric(0.5*crossprod(APY[[i]],PAPY[[j]])) 
        }
    }
    
    AI <- t(AI)
    
    for(i in 1:length(score)){ # this is the row index
        for(j in 1:i){ # this is the column index
            AI[i,j] <- as.numeric(0.5*crossprod(APY[[i]],PAPY[[j]])) 
        }
    }
    
    print(AI); print(score); message(np)
    return(list(AI = AI, score = score))
}	

.calcAIcovMats <- function(PY, covMatList, Sigma.inv, Sigma.inv_X, Xt_Sigma.inv_X.inv){
    m <- length(covMatList)
    AI <- matrix(NA, nrow =  m, ncol = m)
    score <- rep(NA, m)
    
    for (i in 1:m){
        APY <- crossprod(covMatList[[i]],PY)
        PAPY <- crossprod(Sigma.inv, APY) - crossprod(tcrossprod(Xt_Sigma.inv_X.inv, Sigma.inv_X), crossprod(Sigma.inv_X, APY))
        # PAPY <- Sigma.inv %*% crossprod(covMatList[[i]],PY) - tcrossprod(tcrossprod(Sigma.inv_X, Xt_Sigma.inv_X.inv), t(crossprod(covMatList[[i]],PY)) %*% Sigma.inv_X)

        trPA.part1 <- sum( Sigma.inv * covMatList[[i]] )
        trPA.part2 <- sum(diag( (crossprod( Sigma.inv_X, covMatList[[i]]) %*% Sigma.inv_X) %*% Xt_Sigma.inv_X.inv ))
        trPA <-  trPA.part1 - trPA.part2
        
        score[i] <- as.numeric(-0.5*(trPA - crossprod(PY, APY))) # tr(PA) - YPAPY
        AI[i,i] <- as.numeric(0.5*crossprod(APY, PAPY)) # YPAPAPY
        
        # old, less efficient way for calculating score and AI
        # score[i] <- as.numeric(-0.5*(trPA - crossprod(Y, PAPY))) # tr(PA) - YPAPY
        # AI[i,i] <- as.numeric(0.5*crossprod(PY, crossprod(covMatList[[i]],PAPY))) # YPAPAPY
        
        if((i+1) <= m){
            for(j in (i+1):m){
                AI[i,j] <- as.numeric(0.5*crossprod(PY, crossprod(covMatList[[j]],PAPY))) # YPDPAPY
                AI[j,i] <- AI[i,j]
            }
        }
    }
    return(list(AI = AI, score = score))
}	

.calcAIhetvars <- function(PY, group.idx, Sigma.inv, Sigma.inv_X, Xt_Sigma.inv_X.inv){
    g <- length(group.idx)
    n <- length(PY)
    
    if (g == 1){
        PPY <- crossprod(Sigma.inv, PY) - crossprod(tcrossprod(Xt_Sigma.inv_X.inv, Sigma.inv_X), crossprod(Sigma.inv_X, PY))

        trP.part1 <- sum(diag( Sigma.inv ))
        trP.part2 <- sum(diag( crossprod( Sigma.inv_X) %*% Xt_Sigma.inv_X.inv ))
        trP <-  trP.part1 - trP.part2

        score <- as.numeric(-0.5*(trP - crossprod(PY))) # tr(P) - YPIPY
        AI <- as.numeric(0.5*crossprod(PY, PPY)) # YPIPIPY

    } else{	
        AI <- matrix(NA, nrow =  g, ncol = g)
        score <- rep(NA, g)
        
        for (i in 1:g) {
            ### PIPY <- crossprod(P[group.idx[[i]], ], PY[group.idx[[i]]]) 
            ### score[ i] <- -0.5 * (sum(diag(P)[group.idx[[i]]]) - crossprod(PY[group.idx[[i]]]))

            ## covMati <- Diagonal( x=as.numeric( 1:n %in% group.idx[[i]] ) )
            ## IPY <- covMati %*% PY
            ## PIPY <- Sigma.inv %*% IPY - tcrossprod(tcrossprod(Sigma.inv_X, Xt_Sigma.inv_X.inv), t(IPY) %*% Sigma.inv_X)
            ## trPi.part1 <- sum(diag(Sigma.inv)[ group.idx[[i]] ] )
            ## trPi.part2 <- sum(diag( (crossprod( Sigma.inv_X, covMati) %*% Sigma.inv_X) %*% Xt_Sigma.inv_X.inv ))
            
            covMati <- as.numeric( 1:n %in% group.idx[[i]] )
            IPY <- covMati * PY
            PIPY <- crossprod(Sigma.inv, IPY) - crossprod(tcrossprod(Xt_Sigma.inv_X.inv, Sigma.inv_X), crossprod(Sigma.inv_X, IPY))
            
            trPi.part1 <- sum(diag(Sigma.inv)[ group.idx[[i]] ] )
            trPi.part2 <- sum(diag( crossprod(Sigma.inv_X*covMati, Sigma.inv_X) %*% Xt_Sigma.inv_X.inv ))
            trPi <- trPi.part1 - trPi.part2
            
            score[i] <- as.numeric(-0.5*(trPi - crossprod(PY[group.idx[[i]]])))
            AI[i,i] <- as.numeric(0.5*crossprod(PY[group.idx[[i]]], PIPY[group.idx[[i]]]))
            
            if ((i + 1) <= g) {
                for (j in (i + 1):g) {
                    AI[i,j] <- as.numeric(0.5*crossprod(PY[group.idx[[j]]], PIPY[group.idx[[j]]]))
                    AI[j,i] <- AI[i,j]
                }
            }
        }
    }
    return(list(AI = AI, score = score))
}

.calcAIcovMatsResids <- function(PY, covMatList, group.idx, Sigma.inv, Sigma.inv_X, Xt_Sigma.inv_X.inv){
    m <- length(covMatList)
    g <- length(group.idx)
    AI <- matrix(NA, nrow =  m, ncol = g)
    
    for(i in 1:m){
        APY <- crossprod(covMatList[[i]],PY)
        PAPY <- crossprod(Sigma.inv, APY) - crossprod(tcrossprod(Xt_Sigma.inv_X.inv, Sigma.inv_X), crossprod(Sigma.inv_X, APY))
        # PAPY <- Sigma.inv %*% crossprod(covMatList[[i]],PY) - tcrossprod(tcrossprod(Sigma.inv_X, Xt_Sigma.inv_X.inv), t(crossprod(covMatList[[i]],PY)) %*% Sigma.inv_X)   

        ## commenting out these lines for multi-trait analyses for now
##        if(g == 1){
##            AI[i,1] <- as.numeric(0.5*crossprod(PY, PAPY)) # YPIPAPY
##        }else{
##            for(j in 1:g){
##                AI[i,j] <- as.numeric(0.5*crossprod(PY[group.idx[[j]]], PAPY[group.idx[[j]]])) # YP(I_group)PAPY
##            }
##        }
    }

    return(AI)
}
