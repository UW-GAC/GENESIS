pcrelateMakeGRM <- function(pcrelObj, scan.include = NULL, scaleKin = 2){
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

    kinMat <- scaleKin*readex.gdsn(index.gdsn(pcrelObj, "kinship"), sel=list(sample.idx, sample.idx))
    rownames(kinMat) <- scan.include
    colnames(kinMat) <- scan.include

    return(kinMat)
}


