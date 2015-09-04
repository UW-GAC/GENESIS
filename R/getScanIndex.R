getScanIndex <- function(genoData, scan.include){
    # load Scan IDs
    scanIDs <- getScanID(genoData)    

    # samples to include
    if(!is.null(scan.include)){
        if(!all(scan.include %in% scanIDs)){ stop("Not all of the ScanID in scan.include are in genoData") }
    }else{
        scan.include <- scanIDs
    }

    # sample size
    n <- length(scan.include)
    if(n == 0){  stop("None of the samples in scan.include are in genoData")  }

    # get matching index of sample numbers to include
    scan.include.idx <- match(scan.include, scanIDs)
    # get the order of samples in genoData
    scanorder <- order(scan.include.idx)
    # reorder
    scan.include <- scan.include[scanorder]
    scan.include.idx <- scan.include.idx[scanorder]    

    return(list(value = scan.include, index = scan.include.idx, n=n))
}