pcrelateMakeGRM <- function(pcrelObj, scan.include = NULL, scaleKin = 2){
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
    
    # create GRM
    if(class(pcrelObj) == "gds.class"){
        kinMat <- scaleKin*readex.gdsn(index.gdsn(pcrelObj, "kinship"), sel=list(sample.idx, sample.idx))
        id <- readex.gdsn(index.gdsn(pcrelObj, "sample.id"), sel=list(sample.idx))
        rownames(kinMat) <- id
        colnames(kinMat) <- id
    }else if(class(pcrelObj) == "pcrelate"){
        kinMat <- scaleKin*pcrelObj$kinship[sample.idx, sample.idx]
    }
    
    return(kinMat)
}


