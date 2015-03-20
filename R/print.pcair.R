print.pcair <-
function(x,...){
	cat("Call:\n")
	print(x$call)
	cat("\nPCA Method:\n")
	cat(x$method, "\n")
}
