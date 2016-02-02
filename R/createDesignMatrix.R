createDesignMatrix <- function(scanData, outcome, covars, scan.include){
    # check for variable names in scanData
    if(!("scanID" %in% names(scanData))){
        stop("scanID must be in scanData")
    }
    if(!(outcome %in% names(scanData))){
        stop("outcome must be in scanData")
    }

    # set which samples to keep
    scanID <- scanData[,"scanID"]
    # samples to include for the analysis
    if(!is.null(scan.include)){
        if(!all(scan.include %in% scanID)){
            stop("Not all of the scanID in scan.include are in scanData")
        }
        keep <- scanID %in% scan.include
    }else{
        keep <- rep(TRUE, length(scanID))
    }

    # set up model and read in data
    if(!is.null(covars)){
        cvnames <- unique(unlist(strsplit(covars,"[*:]")))
        if(!(all(cvnames %in% names(scanData)))){
            stop("All of the variables in covars must be in scanData")
        }
        dat <- scanData[,c(outcome, cvnames)]
        model.formula <- as.formula(paste(paste(outcome,"~"), paste(covars,collapse="+")))
    }else{
        dat <- scanData[,outcome, drop=FALSE]
        model.formula <- as.formula(paste(paste(outcome,"~"), 1))
    }    
        
    # identify samples with any missing data
    keep <- keep & apply(dat,1,function(x){ all(!is.na(x)) })
    # remove samples with any missing data
    dat <- dat[keep,,drop=FALSE]

    # outcome vector
    Y <- dat[,outcome]
    # create design matrix    
    W <- model.matrix(model.formula, data=dat)
    rownames(W) <- scanID[keep]
    # check for columns of all the same value (except the intercept)
    dropcol <- append(FALSE, apply(W[,-1,drop=FALSE], 2, var) == 0)
    if(sum(dropcol) > 0){
        message("Covariates ",paste(colnames(W)[dropcol], collapse = ", "), " have only 1 value: they have been removed from the model")
        W <- W[,!dropcol,drop=FALSE]
    }
    k <- ncol(W)

    # return output
    return(list(Y = Y, W = W, k = k))
}
