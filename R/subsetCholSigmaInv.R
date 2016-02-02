subsetCholSigmaInv <- function(cholSigmaInv, chol.idx) {
  if(length(chol.idx) > 0){
    # subset cholSigmaInv
    SigmaInv <- tcrossprod(cholSigmaInv)
    for(i in sort(chol.idx, decreasing=TRUE)){
      SigmaInv <- SigmaInv[-i,-i] - tcrossprod(SigmaInv[-i,i])/SigmaInv[i,i]
    }
    cholSigmaInv <- t(chol(SigmaInv))
  }
  
  cholSigmaInv
}
