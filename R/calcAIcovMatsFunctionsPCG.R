.calcAIcovMatsPCG <- function(Y, X, PY, covMatList, Sigma, Sigma.inv_X, Method){
    m <- length(covMatList)
    n <- length(PY)
    AI <- matrix(NA, nrow =  m, ncol = m)
    score <- rep(NA, m)
    
    for (i in 1:m){
        ### PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
        APY <- crossprod(covMatList[[i]],PY)
        PAPY <- .LeftMP(APY, X, Sigma, Sigma.inv_X)
        trPA <-  .LeftMPtrEst(covMatList[[i]], X, Sigma, Sigma.inv_X,Method = Method)
        ###test section
        #A_APYY <- crossprod(covMatList[[i]],(diag(n)-PY %*% t(Y)))
        #score[i] <- as.numeric(-0.5*.LeftMPtrEst(A_APYY,X, Sigma, Sigma.inv_X,Method = Method))
        # print("original var")
        # .theoryvar(covMatList[[i]],X, Sigma, Sigma.inv_X)
        # print("new var")
        #.theoryvar(A_APYY,X, Sigma, Sigma.inv_X)
        #.testtrest(covMatList[[i]], X, Sigma, Sigma.inv_X,Method = Method)
        ### score[i] <- -0.5*(sum(P*covMatList[[i]]) - crossprod(Y, PAPY)) # tr(PA) - YPAPY
        ##test section ends
        score[i] <- as.numeric(-0.5*(trPA - crossprod(Y, PAPY))) # tr(PA) - YPAPY
        AI[i,i] <- as.numeric(0.5*crossprod(PY, crossprod(covMatList[[i]],PAPY))) # YPAPAPY
        if((i+1) <= m){
            for(j in (i+1):m){
                AI[i,j] <- as.numeric(0.5*crossprod(PY, crossprod(covMatList[[j]],PAPY))) # YPDPAPY
                AI[j,i] <- AI[i,j]
            }
        }
    }
    return(list(AI = AI, score = score))
}	


.calcAIcovMatsResidsPCG <- function(Y, X, PY, covMatList, Sigma, Sigma.inv_X,group.idx){
    m <- length(covMatList)
    g <- length(group.idx)
    AI <- matrix(NA, nrow =  m, ncol = g)
    
    for(i in 1:m){
        ### PAPY <- crossprod(P,crossprod(covMatList[[i]],PY))
        APY <- crossprod(covMatList[[i]],PY)
        PAPY <- .LeftMP(APY, X, Sigma, Sigma.inv_X)
        
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


.calcAIhetvarsPCG <- function(Y, X, PY, covMatList, Sigma, Sigma.inv_X,group.idx, Method){
    g <- length(group.idx)
    n <- length(PY)
    
    if (g == 1){
        ### score <-  -0.5*(sum(diag(P)) - crossprod(PY)) 
        ### AI <- 0.5*crossprod(PY,crossprod(P,PY))
        PPY <- .LeftMP(PY, X, Sigma, Sigma.inv_X)
        I <- diag(n)
        trP <-  .LeftMPtrEst(I, X, Sigma, Sigma.inv_X,Method = Method)
        #.theoryvar(I,X, Sigma, Sigma.inv_X)
        #.testtrest(I, X, Sigma, Sigma.inv_X,Method = Method)
        #I_IPYY <- crossprod(I,diag(n) - PY %*% t(Y))
        #score <- -0.5*(.LeftMPtrEst(I_IPYY, X, Sigma, Sigma.inv_X,Method = Method)) #tr(PI) - YPIPY
        
        score <- as.numeric(-0.5*(trP - crossprod(PY)))
        YPIPIPY <- crossprod(PY, PPY)
        AI <- 0.5*as.numeric(YPIPIPY) 
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
            
            I <- as.numeric( 1:n %in% group.idx[[i]] )
            IPY <- I * PY #Not use covMati in the code above -- TIANYU
            PIPY <- .LeftMP(IPY, X, Sigma, Sigma.inv_X)
            #I_IPYY <- crossprod(I,diag(n) - PY %*% t(Y))
            #score[i] <- -0.5*(.LeftMPtrEst(I_IPYY, X, Sigma, Sigma.inv_X,Method = Method)) #tr(PI) - YPIPY
            trPI <- .LeftMPtrEst(I, X, Sigma, Sigma.inv_X,Method = Method)
            #.testtrest(I, X, Sigma, Sigma.inv_X,Method = Method)
            #.theoryvar(I,X, Sigma, Sigma.inv_X)
            score[i] <- -0.5*(trPI - crossprod(PY[group.idx[[i]]])) #tr(PI) - YPIPY
            AI[i,i] <- 0.5 * crossprod(PY[group.idx[[i]]], PIPY[group.idx[[i]]])
            
            if ((i + 1) <= g) {
                for (j in (i + 1):g) {
                    AI[i,j] <- 0.5 * crossprod(PY[group.idx[[j]]], PIPY[group.idx[[j]]]) 
                    AI[j,i] <- AI[i,j]
                }
            }
        }
    }
    return(list(AI = AI, score = score))
}
