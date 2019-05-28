# pchisqsum<- function (x, df, a, remainder="warn", ...)
# {
#     if(!requireNamespace("survey")) stop("package 'survey' must be installed to calculate p-values for SKAT")
#     if(!requireNamespace("CompQuadForm")) stop("package 'CompQuadForm' must be installed to calculate p-values for SKAT")
    
#     ## can happen with randomised trace estimator if most remaining singular values are very small.
#     ##
#     if (any(bad.df <- (df<1))){
#       if (remainder=="warn")
#             warning("Negative/fractional df removed")
#       else if (remainder=="error")
# 	    stop("Negative/fractional df")

#       if(remainder=="missing"){
#             warning("NaN produced from negative/fractional df")
#             return(NaN*x)
# 	}

#       df[bad.df]<-1
#       a[bad.df]<-0
#     }
#     df<-round(df)
#     ##

#     survey::pchisqsum(x, df, a, ...)
# }


pchisqsum_rsvd <- function(x,M,n=100,p=10,q=2, tr2.sample.size=100){
	Q<-srfht2(M,n+p,q=q)
	B<-M%*%Q
	ee<-svd(B,nu=0,nv=0)$d[1:n]^2
	diags <- colSums(M^2)
	tr<-sum(diags)
	if (tr2.sample.size>0){
		tr2<-tracefht(M,k=tr2.sample.size,trace.full=tr)
	} else {
		Ms<-crossprod(M)
		tr2<- sum(Ms^2)
	}	
	tr.small<-tr-sum(ee)
	tr2.small<-tr2-sum(ee^2)
	scale<-tr2.small/tr.small
	nu<-(tr.small^2)/tr2.small
    .pchisqsum(x, c(rep(1,n), nu), c(ee, scale))
}

pchisqsum_ssvd <- function(x,M,n=100,p=10,q=0){
	Q<-srfht(M,n+p,q=q)
	B<-t(Q)%*%M%*%Q
	ee<-eigen(B,symmetric=TRUE,only.values=TRUE)$values[1:n]
	tr<-sum(diag(M))
	tr2<-sum(M^2)
	tr.small<-tr-sum(ee)
	tr2.small<-tr2-sum(ee^2)
	scale<-tr2.small/tr.small
	nu<-(tr.small^2)/tr2.small
    .pchisqsum(x, c(rep(1,n),ceiling(nu)), c(ee, scale))
}

srfht <-
function(A,k,q=0){
	m<-NROW(A)
	mbig<-2^ceiling(log2(m))
	R<-sample(c(-1,1),m,replace=TRUE)
	Astar<-A*R
	idx<-sample(mbig,k)
	AOmega<-fht(Astar)[idx,]/sqrt(k)
	Q<-qr.Q(qr(t(AOmega)))
	for(i in seq_len(q)){
		tildeQ<-qr.Q(qr(crossprod(A,Q)))
		Q<-qr.Q(qr(A%*%tildeQ))
	}
	Q
	}

srfht2 <-
function(A,k,q=0){
	m<-NROW(A)
	mbig<-2^ceiling(log2(m))
	R<-sample(c(-1,1),m,replace=TRUE)
	Astar<-A*R
	idx<-sample(mbig,k)
	AOmega<-fht(Astar)[idx,]/sqrt(k)
	Q<-qr.Q(qr(t(AOmega)))
	for(i in seq_len(q)){
		tildeQ<-qr.Q(qr(A%*%Q))
		Q<-qr.Q(qr(crossprod(A,tildeQ)))
	}
	Q
	}

tracefht <-
function(A,k,trace.full=NULL){
	m<-NROW(A)
	mbig<-2^ceiling(log2(m))
	R<-sample(c(-1,1),m,replace=TRUE)
	Astar<-A*R
	idx<-sample(mbig,k)
	AOmega<-fht(Astar)[idx,]
	AAOmega<-tcrossprod(AOmega,A)
	if (!is.null(trace.full))
		tr<-sum(rowSums(AOmega*AOmega))/k
	trsquared<-sum(rowSums(AAOmega*AAOmega))/k
	if (is.null(trace.full))
		trsquared
	else
		trsquared*(trace.full/tr)^2 
}

fht<-function (A,big=TRUE) 
{
    m <- NROW(A)
    mbig <- 2^ceiling(log2(m))
    Abig <- matrix(0, nrow = mbig, ncol = NCOL(A))
    Abig[1:m, ] <- A
    if (big)
	    .Call("big_mfwht",Abig)
    else 
	    .C("mfwht", Abig, as.integer(mbig), as.integer(NCOL(Abig)))[[1]]
}
