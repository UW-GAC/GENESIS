getVariantInclude <- function(data, variant.include, chromosome){
    # load variant.id
    variantID <- seqGetData(data, "variant.id")

    # variants to include
    if(!is.null(variant.include)){
        # use variants specified variant.include
        if(!all(variant.include %in% variantID)){ stop("Not all of the variant.id in variant.include are in the provided data") }
    }else{
        if(is.null(chromosome)){
            # use all variants
            variant.include <- variantID
        }else{
            # use variantss in specified chromosomes
            variant.include <- variantID[seqGetData(data, "chromosome") %in% chromosome]
        }
    }

    # set a filter
    seqSetFilter(data, variant.id = variant.include, verbose = FALSE)
    # filter monomorphic SNVs
    ref.freq <- alleleFrequency(data, n=0)
    variant.include <- variant.include[ref.freq != 1]

    # number of variants
    n <- length(variant.include)
    if(n == 0){ stop("None of the variants in variant.include are in the provided data") }    

    return(list(value = variant.include, n = n))
}
