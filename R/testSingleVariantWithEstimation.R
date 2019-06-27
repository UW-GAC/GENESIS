calcGtildeWithW <- function(nullmod, G, r=1){
  ##Calculate Gtilde whose square is GWG

  X <- nullmod$model.matrix
  W <- nullmod$W
  
  XWX.inv <- solve(crossprod(X,crossprod(W,X)))
  Gtilde <- G - X %*% (XWX.inv %*% crossprod(X,crossprod(W,G)))
  Gtilde <- r^(1/2)*crossprod(W^(1/2),Gtilde)
  return(Gtilde)
  }

calW <- function(nullmod){
  ###Calculate W matrix

  Y <- nullmod$outcome
  varComp <- nullmod$varComp
  group.idx <- nullmod$group.idx
  vmu <- nullmod$vmu
  
  m <- 1 #number of covariance matrix
  n <- length(Y)
  if (is.null(vmu)){ ## this means the family is "gaussian"
    if (is.null(group.idx)){
      diagV <- rep(varComp[m+1],n)
    } else{
      
      g <- length(group.idx)
      mylevels <- rep(NA, n)
      
      for(i in 1:g){
        mylevels[group.idx[[i]]] <- i # setting up vector of indicators; who is in which group
      }
      diagV <- (varComp[m+1:g])[mylevels]
    }
    
    W <- diag(1/diagV)
  }
  return(W)
}

estr <- function(nullmod, G,
                  rEstSize = 30){
##Estimating the variance ratio r = GPG/GWG

  Y <- nullmod$outcome
  X <- nullmod$model.matrix
  W <- nullmod$W
  
  ncolG <- dim(G)[2]
  
  index <- sample(x = 1:ncolG,size = rEstSize)
  G <- G[,index] ## We choose 30 genes.
  
  Gtilde1 <- calcGtilde(nullmod, G)
  GPG <- colSums(Gtilde1^2)
  
  Gtilde2 <- calcGtildeWithW(nullmod, G)
  GWG <- colSums(Gtilde2^2)
  
  r <- mean(GPG/GWG)
  return(list(GPG,GWG,r))
}
