prepareGenotype <- function(genoData, snp.read.id, scan.read.id, impute.geno){
    # get genotypes for the block
    if(class(genoData) == "GenotypeData"){
        snp.read.idx <- which(getSnpID(genoData) %in% snp.read.id)
        scan.read.idx <- which(getScanID(genoData) %in% scan.read.id)
        geno <- getGenotypeSelection(genoData, snp = snp.read.idx, scan = scan.read.idx, drop = FALSE, transpose = TRUE)
    }else if(class(genoData) == "SeqVarData"){
        seqSetFilter(genoData, variant.id = snp.read.id, sample.id = scan.read.id, verbose = FALSE)
        geno <- altDosage(genoData)
    }
        
    if(impute.geno){
        # impute missing genotype values
        miss.idx <- which(is.na(geno))
        if(length(miss.idx) > 0){
            # get chromosome for selected SNPs
            if(class(genoData) == "GenotypeData"){
                chromChar <- getChromosome(genoData, index = snp.read.idx, char=TRUE)
            }else if(class(genoData) == "SeqVarData"){
                chromChar <- seqGetData(genoData, "chromosome")
            }
            # get sex for selected samples
            if(.hasSex(genoData)){
                if(class(genoData) == "GenotypeData"){
                    sex <- getSex(genoData, index = scan.read.idx)
                }else if(class(genoData) == "SeqVarData"){
                    sex <- sampleData(genoData)$sex
                }
            }else{
                sex <- NULL
            }
            # calculate allele frequencies
            freq <- alleleFreq(geno = geno, chromChar = chromChar, sex = sex)
            miss.snp.idx <- ceiling(miss.idx/length(scan.read.id))
            geno[miss.idx] <- 2*freq[miss.snp.idx]
        }

    }else{
        # check that each sample either has no or all missing genotype data 
        check <- rowSums(is.na(geno))
        if(!all(check %in% c(0, length(snp.read.id)))){
            stop("genoData has sporadic missingness in block size > 1")
        }

        # get rid of samples with completely missing genotype data for this block
        geno <- as.matrix(geno)[(check == 0), , drop=F]
    }

    geno
}
