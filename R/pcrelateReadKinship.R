pcrelateReadKinship <- function(pcrelObj, scan.include = NULL, ibd.probs = TRUE, kin.thresh = NULL){
    # read in sample.id
    sample.id <- read.gdsn(index.gdsn(pcrelObj, "sample.id"))

    # check that the requested samples are in the output
    if(is.null(scan.include)){
        scan.include <- sample.id
    }else{
        if(!all(scan.include %in% sample.id)){
            stop("Some of the samples in scan.include are not in the pcrelObj sample.id")
        }
    }
    # index of requested samples
    sample.idx <- sample.id %in% scan.include
    # number of samples requested
    nsamp <- length(scan.include)

    # indicator matrix for lower triangle
    lowerTriMat <- lower.tri(matrix(nrow=nsamp, ncol=nsamp))

    # set up a data.frame to store results
    out <- data.frame(  ID1 = rep(scan.include[-nsamp], times=((nsamp-1):1)),
                        ID2 = unlist(lapply(2:nsamp,function(x){ scan.include[x:nsamp] })),
                        nsnp = readex.gdsn(index.gdsn(pcrelObj, "nsnp"), sel=list(sample.idx, sample.idx))[lowerTriMat],
                        kin = readex.gdsn(index.gdsn(pcrelObj, "kinship"), sel=list(sample.idx, sample.idx))[lowerTriMat])

    if(ibd.probs){
        # add in the IBD probabilities
        ibd <- readex.gdsn(index.gdsn(pcrelObj, "ibd.probs"), sel=list(sample.idx, sample.idx))
        out$k2 <- ibd[lowerTriMat]
        k0 <- t(ibd)[lowerTriMat]; rm(ibd)
        out$k1 <- 1 - out$k2 - k0
        out$k0 <- k0       
    }


    # keep only those pairs with kinship >= kin.thresh
    if(!is.null(kin.thresh)){
        out <- out[out$kin >= kin.thresh,]
    }

    #if(truncate){
    #    out$kin[out$kin < 0] <- 0
    #    out$kin[out$kin > 1] <- 1
    #    out$k2[out$k2 < 0] <- 0
    #    out$k2[out$k2 > 1] <- 1
    #    out$k1[out$k1 < 0] <- 0
    #    out$k1[out$k1 > 1] <- 1
    #    out$k0[out$k0 < 0] <- 0
    #    out$k0[out$k0 > 1] <- 1
    #}

    return(out)
}

