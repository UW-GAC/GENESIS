getSnpIndex <- function(data, snp.include, chromosome){
    # load snpID
    if(class(data) == "GenotypeData"){
        snpID <- getSnpID(data)
    }else if(class(data) == "SeqVarData"){
        snpID <- seqGetData(data, "variant.id")
    }

    # snps to include
    if(!is.null(snp.include)){
        # use SNPs specified snp.include
        if(!all(snp.include %in% snpID)){ stop("Not all of the snpID in snp.include are in the provided data") }
    }else{
        if(is.null(chromosome)){
            # use all SNPs
            snp.include <- snpID
        }else{
            # use SNPs in specified chromosomes
            if(class(data) == "GenotypeData"){
                snp.include <- snpID[getChromosome(data) %in% chromosome]
            }else if(class(data) == "SeqVarData"){
                snp.include <- snpID[seqGetData(data, "chromosome") %in% chromosome]
            }
        }
    }

    # number of SNPs
    n <- length(snp.include)
    if(n == 0){ stop("None of the SNPs in snp.include are in the provided data") }

    # get matching index of SNP number to include
    snp.include.idx <- match(snp.include, snpID)
    # get the order of SNPs in the data
    snporder <- order(snp.include.idx)
    # reorder
    snp.include <- snp.include[snporder]
    snp.include.idx <- snp.include.idx[snporder]    

    return(list(value = snp.include, index = snp.include.idx, n=n))
}
