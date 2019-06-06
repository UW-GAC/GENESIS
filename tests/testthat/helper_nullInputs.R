.testNullInputs <- function(n=100, binary=FALSE) {
    set.seed(100); x1 <- rnorm(n)
    set.seed(101); x2 <- rbinom(n, size = 1, prob = 0.5)
    X <- cbind(x0=1, x1, x2)

    set.seed(102); sqrt.cor.mat <- matrix(rnorm(n*n, sd = 0.05),n,n, dimnames=list(1:n, 1:n))
    cor.mat <- crossprod(sqrt.cor.mat)

    if (binary) {
	set.seed(104); random.iid <- rnorm(n)
	random <- crossprod(sqrt.cor.mat*0.05, random.iid)
	expit <- function(x){exp(x)/(1+exp(x))} 
	p <- expit(X %*% c(-1, 0.5, 1) + random) 
	set.seed(105); y <- rbinom(n, size = 1, prob = p)
    } else {
        set.seed(106); y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))
    }
    
    group.idx <- list(G1 = c(1:(n/2)), G2 = c((n/2 + 1):n))
    
    return(list(y=y, X=X, cor.mat=cor.mat, group.idx=group.idx))
}

.asDataFrame <- function(dat) {
    df <- data.frame(scanID=1:length(dat$y), y=dat$y, dat$X, group=NA)
    df$group[dat$group.idx[[1]]] <- "G1"
    df$group[dat$group.idx[[2]]] <- "G2"
    df
}

.testGenoMatrix <- function(n=100) {
    set.seed(107); matrix(rbinom(200*n, size = 2, prob = 0.2), nrow = n, ncol = 200)
}

.testNullmod <- function(n=100, MM=FALSE, binary=FALSE) {
    dat <- .testNullInputs(n=n, binary=binary)
    
    if (!binary) {
        family <- "gaussian"
    }  else {
        family <- "binomial"
    }
    
    if (MM) {
        covMatList <- dat$cor.mat
    } else {
        covMatList <- NULL
    }
    
    .fitNullModel(dat$y, dat$X, covMatList = covMatList, group.idx = dat$group.idx, family=family, verbose=FALSE)
}
