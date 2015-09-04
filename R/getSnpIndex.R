getSnpIndex <- function(genoData, snp.include, Xchr){
    # load snpIDs
    snpIDs <- getSnpID(genoData)

    # snps to include
    if(!is.null(snp.include)){
        if(!all(snp.include %in% snpIDs)){ stop("Not all of the SnpID in snp.include are in genoData") }
    }else{
        if(Xchr){
            # all X chromosome SNPs
            snp.include <- snpIDs[getChromosome(genoData) == XchromCode(genoData)]
        }else{
            # all autosomal SNPs
            snp.include <- snpIDs[getChromosome(genoData) <= 22]
        }
    }

    # number of SNPs
    n <- length(snp.include)
    if(n == 0){ stop("None of the SNPs in snp.include are in genoData") }

    # get matching index of SNP number to include
    snp.include.idx <- match(snp.include, snpIDs)
    # get the order of SNPs in the genoData
    snporder <- order(snp.include.idx)
    # reorder
    snp.include <- snp.include[snporder]
    snp.include.idx <- snp.include.idx[snporder]    

    return(list(value = snp.include, index = snp.include.idx, n=n))
}