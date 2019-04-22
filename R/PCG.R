# .pcg <- function(A, b, maxiter=500, tol=1e-5){
#   x0 <- matrix(rep(0, length(b)))
#   Minv <- 1/diag(A)
#   r0 <- b - parallelCorssProduct(A,x0)
#   #r0 <- b - crossprod(A,x0)
#   z0 <- Minv*r0
#   p0 <- matrix(z0)
# 
#   iter <- 0
#   while(iter < maxiter){
#     iter <- iter + 1
# 
#     #Ap0 <- crossprod(A, p0)
#     Ap0 <- parallelCorssProduct(A,p0)
# 
# 
#     alpha <- as.numeric(crossprod(r0, z0)/crossprod(p0, Ap0))
#     x1 <- x0 + alpha*p0
#     r1 <- r0 - alpha*Ap0
# 
#     if(sum(r1^2) < tol){
#       break
#     }else{
#       z1 <- Minv*r1
#       beta <- as.numeric(crossprod(z1, r1)/crossprod(z0, r0))
#       p1 <- z1 + beta*p0
# 
#       # update
#       x0 <- x1
#       r0 <- r1
#       z0 <- z1
#       p0 <- matrix(p1)
#     }
# 
#   }
#   if(iter > maxiter){
#     print("Quit due to iter greater than maxiter, CG may not converge")
#   }
#   # print(paste0("iter from .pcg is",iter))
# 
#   print(paste0("iter from getPCG1ofSigmaAndVector",iter))
#   return(x1)
# }

.pcg <- function(A, b, maxiter=500, tol=1e-5){
  time1 <- Sys.time()
  x0 <- matrix(rep(0, length(b)))
  Minv <- 1/diag(A)
  r0 <- b - crossp(A,x0)
  #r0 <- b - crossprod(A,x0)
  z0 <- Minv*r0
  p0 <- matrix(z0)

  iter <- 0
  while(iter < maxiter){
    iter <- iter + 1

    #Ap0 <- crossprod(A, p0)
    Ap0 <- crossp(A,p0)
    p0Ap0 <- crossprod(p0, Ap0)
    # if(p0Ap0==0){
    #   if(sum(abs(b))==0){
    #     return(rep(0,length(b)))
    #   }
    # }
    alpha <- as.numeric(crossprod(r0, z0)/p0Ap0)
    x1 <- x0 + alpha*p0
    r1 <- r0 - alpha*Ap0

    if(sum(r1^2) < tol){
      break
    }else{
      z1 <- Minv*r1
      beta <- as.numeric(crossprod(z1, r1)/crossprod(z0, r0))
      p1 <- z1 + beta*p0

      # update
      x0 <- x1
      r0 <- r1
      z0 <- z1
      p0 <- matrix(p1)
    }

  }
  if(iter > maxiter){
    print("Quit due to iter greater than maxiter, CG may not converge")
  }
   print(paste0("iter from Rcpp .pcg is",iter))
  time2 <- Sys.time()
  timeonPCG <<- timeonPCG + (time2 - time1)
  return(x1)
}

# 
# .pcg <- function(A, b, maxiter=500, tol=1e-5){
#   x0 <- rep(0, length(b))
#   Minv <- 1/diag(A)
# 
#   r0 <- b - crossprod(A,x0)
#   z0 <- Minv*r0
#   p0 <- z0
# 
#   iter <- 0
#   while(iter < maxiter){
#     iter <- iter + 1
#     Ap0 <- crossprod(A, p0)
#     p0Ap0 <- crossprod(p0, Ap0)
#     if(p0Ap0==0){
#       if(sum(abs(b))==0){
#         return(rep(0,length(b)))
#         }
#       }
#     alpha <- as.numeric(crossprod(r0, z0)/p0Ap0)
#     x1 <- x0 + alpha*p0
#     r1 <- r0 - alpha*Ap0
# 
#     if(sum(r1^2) < tol){
#       break
#     }else{
#       z1 <- Minv*r1
#       beta <- as.numeric(crossprod(z1, r1)/crossprod(z0, r0))
#       p1 <- z1 + beta*p0
# 
#       # update
#       x0 <- x1
#       r0 <- r1
#       z0 <- z1
#       p0 <- p1
#     }
# 
#   }
#   if(iter > maxiter){
#     print("Quit due to iter greater than maxiter, CG may not converge")
#   }
#   # print(paste0("iter from .pcg is",iter))
# 
#   return(x1)
# }

.pcgM <- function(A, M, maxiter=1e5, tol=1e-6){
  # print(".pcgM is running")
  # print(paste0("the dimension of matrix M is ",dim(M)[2]))
  # if(dim(M)[2]>9){
  # print(M)}
  if(!is.null(dim(M))){
    AinvM <- matrix(0, ncol=dim(M)[2], nrow=dim(M)[1])
    for(i in 1:dim(M)[2]){
      AinvM[,i] <- .pcg(A, M[,i], maxiter, tol)
    }
    # print(".pcgM finishes")
    return(AinvM)
  }
  else{
    # print(".pcgM is finishes")
    return(.pcg(A,M))
  }
}

.LeftMPtrEst <- function(Y, X, Sigma, Sigma.inv_X, Method, n=30){
  ##ESTIMATE THE TRACE OF PY.
  # print("Trace estimation begins...")
  if(Method == "EST"){
  Ytrh <- trh <- matrix(0, ncol = 1, nrow = n)
  for (i in 1:n){
    if(!is.null(dim(Y))){
      e <- 2*rbinom(dim(Y)[1],1,0.5)-1
    }
    else{
     e <- 2*rbinom(length(Y),1,0.5)-1
    }
    Ye <- Y %*% e
    PYe <- .LeftMP(Ye, X, Sigma, Sigma.inv_X)
    trh[i] <- crossprod(e,PYe)
    Ytrh[i] <- crossprod(e,Ye)
  }
  trest <- mean(trh)
  return(trest)
  if(var(Ytrh)<10^(-4)) return(trest)
  trYest <- mean(Ytrh)
  trY <- sum(diag(Y))
  a <- cor(Ytrh,trh)
  trest.stab <- trest - a*(trYest - trY)
  # trPY <- sum(diag(.LeftMP(Y, X, Sigma, Sigma.inv_X)))
  # old.diff <<- c(old.diff,trest - trPY)
  # new.diff <<- c(new.diff,trest.stab - trPY)
  return(trest.stab)
  }
  if(Method == "EXT"){
    #print(sum(diag(.LeftMP(Y, X, Sigma, Sigma.inv_X))))
    # print("Trace estimation ends")
    return(sum(diag(.LeftMP(Y, X, Sigma, Sigma.inv_X))))
    }
}

.LeftMP <- function(Y,X, Sigma, Sigma.inv_X){
  #Y LEFT MULTIPLIED BY P MATRIX
  # print("now calculate sigma-1Y")
  # print(Sys.time())
  Sigma.inv_Y <- .pcgM(Sigma, Y)
  # print("now finish sigma-1Y")
  # print(Sys.time())
  Mid <- crossprod(X,Sigma.inv_X)
  Right <- crossprod(X,Sigma.inv_Y)
  # print("now calculate mid,right")
  # print(Sys.time())
   MidRight <- .pcgM(Mid, Right)
  # print("now finish mid,right")
  # print(Sys.time())
  PY <- Sigma.inv_Y - Sigma.inv_X %*% MidRight
  return(PY)
}

.computeSigma <- function(varComp, covMatList, group.idx = NULL, vmu = NULL, gmuinv = NULL){
  m <- length(covMatList)
  n <- nrow(covMatList[[1]])
  
  ###    Sigma <- Vre <- Reduce("+", mapply("*", covMatList, varComp[1:m], SIMPLIFY=FALSE))
  Vre <- Reduce("+", mapply("*", covMatList, varComp[1:m], SIMPLIFY=FALSE))
  ##tau psi
  #m=1 for only kinship
  ##\sum_i \sigma_i^2M_i
  
  if (is.null(vmu)){ ## this means the family is "gaussian"
    if (is.null(group.idx)){
      #Sigma <- Diagonal(x=rep(varComp[m+1],n)) + Vre
      diagV <- rep(varComp[m+1],n)
      ## don't allow
    } else{
      g <- length(group.idx)
      ###        diagV <- rep(0,nrow(covMatList[[1]]))
      ###        for(i in 1:g){
      ###            diagV[group.idx[[i]]] <- varComp[m+i]
      ###        }
      ###        diag(Sigma) <- diag(Sigma) + diagV
      
      mylevels <- rep(NA, n)
      for(i in 1:g){
        mylevels[group.idx[[i]]] <- i # setting up vector of indicators; who is in which group
      }
      #Sigma <- Diagonal(x=(varComp[m+1:g])[mylevels] ) + Vre
      diagV <- (varComp[m+1:g])[mylevels]
      ##allow residual variance to be diff.
      ##W^{-1}
    }
    
    ### if non-gaussian family:
  } else {
    ###        Sigma <- Sigma + diag(as.vector(vmu)/as.vector(gmuinv)^2)
    #Sigma <- Diagonal(x=as.vector(vmu)/as.vector(gmuinv)^2) + Vre
    diagV <- as.vector(vmu)/as.vector(gmuinv)^2
  }
  
  # add diagonal elements
  Sigma <- Vre
  diag(Sigma) <- diag(Sigma) + diagV

  return(list(Sigma = Sigma,Vre = Vre))
  
}

.calREML <- function(Y, X, Sigma, Sigma.inv_X, PY){
  n <- length(Y)
  k <- ncol(X)
  cR <- -0.5*(n-k)*log(2*pi)
  det1 <- -0.5*log(det(Sigma))
  det2 <- -0.5*log(det(crossprod(X,Sigma.inv_X)))
  REML <- as.numeric(cR + det1 + det2 -0.5*crossprod(Y, PY))
  print(REML)
  return(REML)
  }
.testtrest <- function(Y,X, Sigma, Sigma.inv_X, Method, N = 100){
  est <- matrix(0,N)
  for(j in 1:N){
    est[j] <- .LeftMPtrEst(Y, X, Sigma, Sigma.inv_X, Method)
  }
  # print(summary(est))
  # print(sum(diag(.LeftMP(Y,X, Sigma, Sigma.inv_X))))
  print("CV IS")
  print(sd(est)/mean(est))
  }

.theoryvar <- function(Y,X, Sigma, Sigma.inv_X){
  M <- .LeftMP(Y,X, Sigma, Sigma.inv_X)
  varr <- 2*(norm(M,type = "F")^2 - sum(diag(M)^2))
  #print("theory var is")
  print(varr)
  }