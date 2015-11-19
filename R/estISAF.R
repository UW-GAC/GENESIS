estISAF <- function(geno, pcMat, pcProd, transpose = FALSE){
    # rather than compute new hat matrix for each SNP, impute missing genotype values to sample mean
    
    # which genotype values are missing
    miss.idx <- which(is.na(geno))
    # if there are missing genotype values
    if(length(miss.idx) > 0){
        # calculate MAF
        pA <- 0.5*rowMeans(geno, na.rm = TRUE)
        # number of snps in the block
        nsnp.block <- nrow(geno)
        # index of which snps have the missing values
        snp.idx <- miss.idx %% nsnp.block; snp.idx[snp.idx==0] <- nsnp.block
        # replace missing genotypes with twice the sample allele frequency
        geno[miss.idx] <- 2*pA[snp.idx]
    }

    # matrix of individual specific allele frequencies
    if(transpose){
        muhat <- 0.5*tcrossprod(pcMat, tcrossprod(geno, pcProd) )
    }else{
        muhat <- 0.5*tcrossprod( tcrossprod(geno, pcProd), pcMat)
    }
    return(muhat)
}
