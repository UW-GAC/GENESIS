## idx.exclude are indices of individuals that should be excluded (e.g. because of missing genotypes)
nullModelSubset <- function(nullmod, idx.exclude){
    nullmod$fit <- nullmod$fix[-idx.exclude, ]
    nullmod$model.matrix <- nullmod$model.matrix[-idx.exclude,]

    if (nullmod$family$mixedmodel){  ## n by n cholSigmaInv {
        nullmod$cholSigmaInv <- subsetCholSigmaInv(nullmod$cholSigmaInv, idx.exclude)
    }

    if (!nullmod$family$mixedmodel & (nullmod$family$family == "gaussian")){  ## a diagonal or scalar cholSigmaInv

        if (nullmod$hetResid)	{  ## cholSigmaInv is diagonal
            nullmod$cholSigmaInv <- nullmod$cholSigmaInv[-idx.exclude, -idx.exclude]
        }
    }
}


# this is a fancy way of getting the inverse of the subset without having to get the original matrix
# cholesky decomposition of sigma inverse (inverse phenotype covariance matrix)
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
