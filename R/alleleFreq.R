alleleFreq <- function(geno, chromChar, sex, sample.use){
    if(missing(sample.use)){
        sample.use <- rep(TRUE, nrow(geno))
    }
    
    # type of chromosome for each SNP
    xchr <- chromChar == "X"
    ychr <- chromChar == "Y"
    auto <- !xchr & !ychr

    # allele frequency vector
    freq <- rep(NA, ncol(geno))

    # autosomes
    freq[auto] <- 0.5*colMeans(geno[sample.use, auto, drop=FALSE], na.rm=TRUE)

    # xchr
    if(sum(xchr) > 0){
        female <- sex == "F" & sample.use
        male <- sex == "M" & sample.use
        F.count <- colSums(geno[female, xchr, drop = FALSE], na.rm = TRUE)
        F.nsamp <- colSums(!is.na(geno[female, xchr, drop = FALSE]))
        M.count <- 0.5*colSums(geno[male, xchr, drop = FALSE], na.rm = TRUE)
        M.nsamp <- colSums(!is.na(geno[male, xchr, drop = FALSE]))
        freq[xchr] <- (F.count + M.count)/(2*F.nsamp + M.nsamp)
    }
    
    # ychr
    if(sum(ychr) > 0){
        male <- sex == "M" & sample.use
        freq[ychr] <- 0.5*colMeans(geno[male, ychr, drop = FALSE], na.rm=TRUE)
    }

    freq
}
