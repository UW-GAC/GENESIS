pcrelateReadKinship <- function(pcrelObj, scan.include = NULL, ibd.probs = TRUE, kin.thresh = NULL){
    # read in sample.id
    if(class(pcrelObj) == "gds.class"){
        sample.id <- read.gdsn(index.gdsn(pcrelObj, "sample.id"))
    }else if(class(pcrelObj) == "pcrelate"){
        sample.id <- pcrelObj$sample.id
    }

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
    out <- data.frame(ID1 = rep(scan.include[-nsamp], times=((nsamp-1):1)),
                      ID2 = unlist(lapply(2:nsamp,function(x){ scan.include[x:nsamp] })))
    
    if(class(pcrelObj) == "gds.class"){
        out$nsnp <- readex.gdsn(index.gdsn(pcrelObj, "nsnp"), sel=list(sample.idx, sample.idx))[lowerTriMat]
        out$kin <- readex.gdsn(index.gdsn(pcrelObj, "kinship"), sel=list(sample.idx, sample.idx))[lowerTriMat]
    }else if(class(pcrelObj) == "pcrelate"){
        out$nsnp <- pcrelObj$nsnp[sample.idx, sample.idx][lowerTriMat]
        out$kin = pcrelObj$kinship[sample.idx, sample.idx][lowerTriMat]
    }

    if(ibd.probs){
        # add in the IBD probabilities
        if(class(pcrelObj) == "gds.class"){
            ibd <- readex.gdsn(index.gdsn(pcrelObj, "ibd.probs"), sel=list(sample.idx, sample.idx))
        }else if(class(pcrelObj) == "pcrelate"){
            ibd <- pcrelObj$ibd.probs[sample.idx, sample.idx]
        }
        out$k0 <- ibd[lowerTriMat]
        k2 <- t(ibd)[lowerTriMat]; rm(ibd)
        out$k1 <- 1 - out$k0 - k2
        out$k2 <- k2
    }


    # keep only those pairs with kinship >= kin.thresh
    if(!is.null(kin.thresh)){
        out <- out[out$kin >= kin.thresh,]
    }

    return(out)
}

