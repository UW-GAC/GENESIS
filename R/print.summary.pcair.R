print.summary.pcair <-
function(x,...){
	cat("Call:\n")
	print(x$call)
	cat("\nPCA Method:", x$method, "\n")
	if(x$method == "PC-AiR"){
	    cat("\nSample Size:", x$nsamp, "\n")
	    cat("Unrelated Set:", length(x$unrels), "Samples \n")
	    cat("Related Set:", length(x$rels), "Samples \n")
	    cat("\nKinship Threshold:", x$kin.thresh, "\n")
	    cat("Divergence Threshold:", -abs(x$div.thresh), "\n")
	}else{
	    cat("\nSample Size:", x$nsamp, "\n")
	}
	cat("\nPrincipal Components Returned:", dim(x$vectors)[2], "\n")
	cat("Eigenvalues:", round(x$values[1:min(10,length(x$values))],3), "...\n")
	cat("\nMAF Filter:", x$MAF, "\n")
	cat("SNPs Used:", x$nsnps, "\n")	
}
