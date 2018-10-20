setGeneric("makeSparseMatrix", function(x, ...) standardGeneric("makeSparseMatrix"))

setMethod("makeSparseMatrix",
          "matrix",
          function(x, thresh = NULL, sample.include = NULL, diag.value = NULL, verbose = TRUE){
          	# melt to a data.table
          	x <- meltMatrix(x = x, drop.lower = TRUE, drop.diag = ifelse(is.null(diag.value), FALSE, TRUE))
          	.makeSparseMatrix(x = x, thresh = thresh, sample.include = sample.include, diag.value = diag.value, verbose = verbose)
          })

setMethod("makeSparseMatrix",
          "Matrix",
          function(x, thresh = NULL, sample.include = NULL, diag.value = NULL, verbose = TRUE){
          	# melt to a data.table
          	x <- meltMatrix(x = x, drop.lower = TRUE, drop.diag = ifelse(is.null(diag.value), FALSE, TRUE))
          	.makeSparseMatrix(x = x, thresh = thresh, sample.include = sample.include, diag.value = diag.value, verbose = verbose)
          })

setMethod("makeSparseMatrix",
          "data.frame",
          function(x, thresh = NULL, sample.include = NULL, diag.value = NULL, verbose = TRUE){
          	.makeSparseMatrix(x = as.data.table(x), thresh = thresh, sample.include = sample.include, diag.value = diag.value, verbose = verbose)
          })

setMethod("makeSparseMatrix",
          "data.table",
          function(x, thresh = NULL, sample.include = NULL, diag.value = NULL, verbose = TRUE){
            .makeSparseMatrix(x = x, thresh = thresh, sample.include = sample.include, diag.value = diag.value, verbose = verbose)
          })



.makeSparseMatrix <- function(x, thresh = NULL, sample.include = NULL, diag.value = NULL, verbose = TRUE){

    # keep R CMD check from warning about undefined global variables
    ID1 <- ID2 <- value <- NULL
    `.` <- function(...) NULL

	# check sample.include
	if(!is.null(sample.include)){
        # subset to samples in sample.include
        if(verbose) message('Using ', length(sample.include), ' samples in sample.include')
        x <- x[ID1 %in% sample.include & ID2 %in% sample.include]
    }else{
        # get list of all samples in the data
        sample.include <- sort(unique(c(x$ID1, x$ID2)))
        if(verbose) message("Using ", length(sample.include), " samples provided")
    }

    # check for diag values
    if(is.null(diag.value)){
    	if(!all.equal(sort(x[ID1 == ID2, ID1]), sort(sample.include))) stop('When `diag.value` is NULL, diagonal values must be provided for all samples')
    }
    
    # create graph of relatives
    if(verbose) message("Identifying clusters of relatives...")
    if(is.null(thresh)){
        g <- igraph::graph_from_data_frame(x[ID1 != ID2])
    }else{
        g <- igraph::graph_from_data_frame(x[ID1 != ID2 & value > thresh])
    }
    # extract cluster membership
    clu <- igraph::components(g)
    mem <- clu$membership
    if(verbose) message("    ", length(mem), " relatives in ", clu$no, " clusters")
    
    blocks <- list()
    block.id <- list()
    if(clu$no > 0){
        if(verbose) message("Creating block matrices for clusters...")
        for(i in 1:clu$no){
            # samples in the cluster
            ids <- names(mem[mem == i])
            # create a table for all pairs in the cluster
            allpairs <- as.data.table(expand.grid(ID2 = ids, ID1 = ids))
            
            # merge
            sub <- x[ID1 %in% ids & ID2 %in% ids][allpairs, on = c("ID1", "ID2")]
            # set pairs without values to 0
            sub[is.na(value), value := 0]
            
            # cast to a matrix
            submat <- reshape2::acast(data = sub, formula = ID1 ~ ID2, value.var = "value")
            # put the values on both sides of the diagonal
            submat <- submat + t(submat)
            # fix the diagonal
            if(is.null(diag.value)){
                diag(submat) <- 0.5*diag(submat)
            }else{
                diag(submat) <- diag.value
            }
            
            # store in the list
            blocks[[i]] <- submat
            block.id[[i]] <- rownames(submat)
        }
    }else{
        if(verbose) message("No clusters identified")
    }
        
    # add in identity matrix of unrelated samples
    unrel.id <- setdiff(sample.include, names(mem))
    if(verbose) message(length(unrel.id), " samples with no relatives included...")

    if(is.null(diag.value)){
        # data for the diagonal
        ddat <- x[ID1 == ID2, .(ID1, value)][data.table(ID1 = unrel.id), on = 'ID1']
        blocks[[clu$no + 1]] <- Diagonal(n = nrow(ddat), x = ddat$value)
        block.id[[clu$no + 1]] <- ddat$ID1
    }else{
        blocks[[clu$no + 1]] <- Diagonal(n = length(unrel.id), x = diag.value)
        block.id[[clu$no + 1]] <- unrel.id
    }

    # create block diagonal matrix
    if(verbose) message("Putting all samples together into one block diagonal matrix")
    mat_sparse <- bdiag(blocks)
    # ids of samples
    mat.id <- unlist(block.id)
    rownames(mat_sparse) <- mat.id
    colnames(mat_sparse) <- mat.id

    # set the matrix to symmetric to save memory
    mat_sparse <- as(mat_sparse, "symmetricMatrix")
    return(mat_sparse)
}



setGeneric("meltMatrix", function(x, ...) standardGeneric("meltMatrix"))

setMethod("meltMatrix",
          "matrix",
          function(x, drop.lower = FALSE, drop.diag = FALSE){
            if(drop.lower){
                x[lower.tri(x, diag = drop.diag)] <- NA
            }
            x <- as.data.table(reshape2::melt(x, varnames = c('ID1', 'ID2'), na.rm = TRUE, as.is = TRUE))
            setkeyv(x, c('ID1', 'ID2'))
          })

setMethod("meltMatrix",
          "Matrix",
          function(x, drop.lower = FALSE, drop.diag = FALSE){
            if(drop.lower){
                x[lower.tri(x, diag = drop.diag)] <- NA
            }

            # this modeled after reshape2:::melt.matrix
            labels <- as.data.table(expand.grid(dimnames(x), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
            setnames(labels, c('Var1', 'Var2'), c('ID1', 'ID2'))

            missing <- is.na(as.vector(x))
            x <- cbind(labels[!missing,], data.table(value = x[!missing])) 
            setkeyv(x, c('ID1', 'ID2'))
          })


pcrelateMakeGRM <- function(pcrelobj, sample.include = NULL, thresh = NULL, scaleKin = 2, verbose = TRUE){
    if(!requireNamespace("igraph")) stop("package 'igraph' must be installed to use this function")

    # get the diagonals
    x <- pcrelobj$kinSelf[, .(ID, f)]
    setnames(x, 'ID', 'ID1')
    x[, ID2 := ID1]
    x[, kin := 0.5*(1 + f)][, f := NULL]

    # append the off-diagonal
    x <- rbind(x, pcrelobj$kinBtwn[,.(ID1, ID2, kin)])

    # scale the values
    x[, kin := scaleKin*kin]
    
    # set kin name to value
    setnames(x, 'kin', 'value')
    # call makeSparseMatrix
    makeSparseMatrix(x = x, thresh = thresh, sample.include = sample.include, diag.value = NULL, verbose = verbose)
}


kingToMatrix <- function(file.king, sample.include = NULL, thresh = NULL, verbose = TRUE){
    if(!requireNamespace("igraph")) stop("package 'igraph' must be installed to use this function")

    # keep R CMD check from warning about undefined global variables
    ID1 <- ID2 <- PropIBD <- Kinship <- NULL
    `.` <- function(...) NULL
    
    # read in the king results and subset columns
    if(verbose) message('Reading in KING output...')
    
    # if multiple input files
    if(length(file.king) > 1){
        king <- lapply(file.king, fread)
        # pull out column names in common
        cnames <- Reduce(intersect, lapply(king, colnames))
        # subset and rbind
        king <- do.call(rbind, lapply(king, function(x){ x[, colnames(x) %in% cnames, with = FALSE]} ))
        
    # one input file
    }else{
        king <- fread(file.king)
    }
    
    # subset to needed columns
    if('PropIBD' %in% colnames(king)){
        if(verbose) message('Inferred to be KING --ibdseg output')
        king <- king[, .(ID1, ID2, PropIBD)]
        king <- king[, value := 0.5*PropIBD][, PropIBD := NULL]
    }else if('Kinship' %in% colnames(king)){
        if(verbose) message('Inferred to be KING --robust output')
        king <- king[, .(ID1, ID2, Kinship)]
        setnames(king, 'Kinship', 'value')
    }else{
        stop('All files in file.king must contain a column called `PropIBD` (KING --ibdseg output) or `Kinship` (KING --robust output)')
    }
        
    # check for duplicate pairs
    setkeyv(king, c('ID1', 'ID2'))
    if(any(duplicated(king))){
        stop('Some sample pairs are provided multiple times in file.king; please only provide one value per sample pair')
    }

    # call makeSparseMatrix
    makeSparseMatrix(x = king, thresh = thresh, sample.include = sample.include, diag.value = 0.5, verbose = verbose)
}
