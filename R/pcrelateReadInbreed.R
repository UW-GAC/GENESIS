pcrelateReadInbreed <- function(pcrelObj, scan.include = NULL, f.thresh = NULL){
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

    # set up a data.frame to store results
    out <- data.frame(  ID = scan.include,
                        nsnp = diag(readex.gdsn(index.gdsn(pcrelObj, "nsnp"), sel = list(sample.idx, sample.idx))),
                        f = 2*diag(readex.gdsn(index.gdsn(pcrelObj, "kinship"), sel=list(sample.idx, sample.idx)))-1)

    # keep only those pairs with inbreeding >= f.thresh
    if(!is.null(f.thresh)){
        out <- out[out$f >= f.thresh,]
    }

    return(out)
}


