prepareGenotype <- function(genoData, snp.read.idx, scan.read.idx, impute.geno){
    # get genotypes for the block
    geno <- getGenotypeSelection(genoData, snp = snp.read.idx, scan = scan.read.idx, drop = FALSE, transpose = TRUE)
    
    if(impute.geno){
        # impute missing genotype values
        miss.idx <- which(is.na(geno))
        if(length(miss.idx) > 0){
            # get chromosome for selected SNPs
            chromChar <- getChromosome(genoData, index = snp.read.idx, char=TRUE)
            # get sex for selected samples
            if(hasSex(genoData)){
                sex <- getSex(genoData, index = scan.read.idx)
            }else{
                sex <- NULL
            }
            # calculate allele frequencies
            freq <- alleleFreq(geno = geno, chromChar = chromChar, sex = sex)
            miss.snp.idx <- ceiling(miss.idx/length(scan.read.idx))
            geno[miss.idx] <- 2*freq[miss.snp.idx]
        }

    }else{
        # check that each sample either has no or all missing genotype data 
        check <- rowSums(is.na(geno))
        if(!all(check %in% c(0, length(snp.read.idx)))){
            stop("genoData has sporadic missingness in block size > 1")
        }

        # get rid of samples with completely missing genotype data for this block
        geno <- as.matrix(geno)[(check == 0), , drop=F]
    }

    geno
}
