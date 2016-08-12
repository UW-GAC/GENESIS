getScanIndex <- function(data, scan.include){    
    # load Scan IDs
    if(class(data) == "ScanAnnotationDataFrame" | class(data) == "GenotypeData"){
        scanID <- getScanID(data)

    }else if(class(data) == "SeqVarData"){
        scanID <- seqGetData(data, "sample.id")
    
    }else if(class(data) == "data.frame"){
        if(!("scanID" %in% names(data))){
            stop("scanID must be in provided data")
        }
        scanID <- data[,"scanID"]
    }
    
    # samples to include
    if(!is.null(scan.include)){
        if(!all(scan.include %in% scanID)){ 
            stop("Not all of the scanID in scan.include are in the provided data") 
        }
    }else{
        scan.include <- scanID
    }

    # sample size
    n <- length(scan.include)
    if(n == 0){  stop("None of the samples in scan.include are in the provided data")  }

    # get matching index of sample numbers to include
    scan.include.idx <- match(scan.include, scanID)
    # get the order of samples in genoData
    scanorder <- order(scan.include.idx)
    # reorder
    scan.include <- scan.include[scanorder]
    scan.include.idx <- scan.include.idx[scanorder]    

    return(list(value = scan.include, index = scan.include.idx, n=n))
}
