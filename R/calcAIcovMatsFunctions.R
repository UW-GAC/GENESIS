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
            AI[i,i] <- as.numeric(0.5*crossprod(PY[group.idx[[i]]], PIPY[group.idx[[i]]])
            
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
        
        if(g == 1){
            AI[i,1] <- as.numeric(0.5*crossprod(PY, PAPY)) # YPIPAPY
        }else{
            for(j in 1:g){
                AI[i,j] <- as.numeric(0.5*crossprod(PY[group.idx[[j]]], PAPY[group.idx[[j]]])) # YP(I_group)PAPY
            }
        }
    }

    return(AI)
}