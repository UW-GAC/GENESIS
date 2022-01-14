setGeneric("makeSparseMatrix", function(x, ...) standardGeneric("makeSparseMatrix"))

setMethod("makeSparseMatrix",
          "matrix",
          function(x, thresh = 2^(-11/2), sample.include = NULL, diag.value = NULL, verbose = TRUE){
            .makeSparseMatrix_matrix(x = x, thresh = thresh, sample.include = sample.include, diag.value = diag.value, verbose = verbose)
          })

setMethod("makeSparseMatrix",
          "Matrix",
          function(x, thresh = 2^(-11/2), sample.include = NULL, diag.value = NULL, verbose = TRUE){
          	.makeSparseMatrix_matrix(x = x, thresh = thresh, sample.include = sample.include, diag.value = diag.value, verbose = verbose)
          })

setMethod("makeSparseMatrix",
          "data.frame",
          function(x, thresh = NULL, sample.include = NULL, diag.value = NULL, verbose = TRUE){
          	.makeSparseMatrix_df(x = as.data.table(x), thresh = thresh, sample.include = sample.include, diag.value = diag.value, verbose = verbose)
          })

setMethod("makeSparseMatrix",
          "data.table",
          function(x, thresh = NULL, sample.include = NULL, diag.value = NULL, verbose = TRUE){
            .makeSparseMatrix_df(x = x, thresh = thresh, sample.include = sample.include, diag.value = diag.value, verbose = verbose)
          })


.makeSparseMatrix_matrix <- function(x, thresh = 2^(-11/2), sample.include = NULL, diag.value = NULL, verbose = TRUE){

    # keep R CMD check from warning about undefined global variables
    ID1 <- ID2 <- NULL

    # check that thresh is not NULL
    if(is.null(thresh)){
        stop("thresh must be specified when x is a matrix object")
    }

    # check for names
    if(is.null(colnames(x))) {
        if(is.null(rownames(x))) {
            stop("dimnames must be provided for matrix x to track sample order in the output matrix.")
        } else {
            colnames(x) <- rownames(x)
        }
    }
    if(is.null(rownames(x))) {
        rownames(x) <- colnames(x)
    }

    # check that matrix is square
    if(!all(rownames(x) == colnames(x))){
        stop("x must be square when it is a matrix object")
    }

    # check sample.include
    if(!is.null(sample.include)){
        # subset to samples in sample.include
        if(verbose) message('Using ', length(sample.include), ' samples in sample.include')
        x <- x[rownames(x) %in% sample.include, colnames(x) %in% sample.include]
    }else{
        # get list of all samples in the data
        sample.include <- sort(rownames(x))
        if(verbose) message("Using ", length(sample.include), " samples provided")
    }

    # check for diag values
    if(is.null(diag.value)){
        if(any(is.na(diag(x)))) stop('When `diag.value` is NULL, diagonal values must be provided for all samples')
    }

    # get the table of all related pairs
    rel <- .apply(x, MARGIN = 1, FUN = function(v){ names(v)[v > thresh & !is.na(v)] })
    rel <- lapply(seq_along(rel), function(i) { data.table('ID1' = names(rel)[[i]], 'ID2' = rel[[i]]) })
    rel <- rbindlist(rel)
    rel <- rel[,`:=`(ID1 = as.character(ID1), ID2 = as.character(ID2))]
    setkeyv(rel, c('ID1', 'ID2'))

    # create graph of relatives
    if(verbose) message("Identifying clusters of relatives...")
    g <- igraph::graph_from_data_frame(rel[ID1 != ID2])
    # extract cluster membership
    clu <- igraph::components(g)
    mem <- clu$membership

    blocks <- list()
    block.id <- list()
    if(clu$no > 0){
        if(verbose) message("    ", length(mem), " relatives in ", clu$no, " clusters; largest cluster = ", max(clu$csize))
        if(verbose) message("Creating block matrices for clusters...")
        for(i in 1:clu$no){
            # samples in the cluster
            ids <- names(mem[mem == i])

            # extract the matrix for all pairs in the cluster
            submat <- x[rownames(x) %in% ids, colnames(x) %in% ids]

            # fix the diagonal
            if(!is.null(diag.value)){
                diag(submat) <- diag.value
            }

            # store in the list
            blocks[[i]] <- submat
            block.id[[i]] <- rownames(submat)
        }
    }else{
        if(verbose) message("    No clusters identified")
    }

    # add in identity matrix of unrelated samples
    unrel.id <- setdiff(sample.include, names(mem))
    if(verbose) message(length(unrel.id), " samples with no relatives included")

    if(length(unrel.id) > 0) {
        if(is.null(diag.value)){
            # data for the diagonal
            ddat <- diag(x)[rownames(x) %in% unrel.id]
            blocks[[clu$no + 1]] <- Diagonal(n = length(ddat), x = ddat)
            block.id[[clu$no + 1]] <- names(ddat)
        }else{
            blocks[[clu$no + 1]] <- Diagonal(n = length(unrel.id), x = diag.value)
            block.id[[clu$no + 1]] <- unrel.id
        }
    }

    # create block diagonal matrix
    if(length(blocks) > 1){
        if(verbose) message("Putting all samples together into one block diagonal matrix")
        mat_sparse <- bdiag(blocks)
    }else{
        mat_sparse <- Matrix(blocks[[1]])
    }

    # ids of samples
    mat.id <- unlist(block.id)
    rownames(mat_sparse) <- mat.id
    colnames(mat_sparse) <- mat.id

    # set the matrix to symmetric to save memory
    mat_sparse <- as(mat_sparse, "symmetricMatrix")
    return(mat_sparse)
}

.makeSparseMatrix_df <- function(x, thresh = NULL, sample.include = NULL, diag.value = NULL, verbose = TRUE){

    # keep R CMD check from warning about undefined global variables
    ID1 <- ID2 <- value <- NULL

    # make sure IDs are character
    if (!is.character(x$ID1)) x$ID1 <- as.character(x$ID1)
    if (!is.character(x$ID2)) x$ID2 <- as.character(x$ID2)

    # check sample.include
    if(!is.null(sample.include)){
        # subset to samples in sample.include
        if(verbose) message('Using ', length(sample.include), ' samples in sample.include')
        x <- x[ID1 %in% sample.include & ID2 %in% sample.include]
    }else{
        # get list of all samples in the data
        sample.include <- unique(rbind(setnames(unique(x[,'ID1']), 'ID1', 'ID'), setnames(unique(x[,'ID2']), 'ID2', 'ID'))[,'ID'])
        setkey(sample.include, 'ID')
        sample.include <- sample.include[['ID']]
        # sample.include <- sort(unique(c(x$ID1, x$ID2)))
        if(verbose) message("Using ", length(sample.include), " samples provided")
    }

    # check for diag values
    if(is.null(diag.value)){
    	if(!setequal(x[ID1 == ID2, ID1], sample.include)) stop('When `diag.value` is NULL, diagonal values must be provided for all samples')
    }

    # create graph of relatives
    if(verbose) message("Identifying clusters of relatives...")
    if(is.null(thresh)){
        g <- igraph::graph_from_data_frame(x[ID1 != ID2 & value != 0])
    }else{
        g <- igraph::graph_from_data_frame(x[ID1 != ID2 & value > thresh])
    }
    # extract cluster membership
    clu <- igraph::components(g)
    mem <- clu$membership

    blocks <- list()
    block.id <- list()
    if(clu$no > 0){
        if(verbose) message("    ", length(mem), " relatives in ", clu$no, " clusters; largest cluster = ", max(clu$csize))
        if(verbose) message("Creating block matrices for clusters...")
        for(i in 1:clu$no){
            # samples in the cluster
            ids <- names(mem[mem == i])
            # create a table for all pairs in the cluster
            allpairs <- as.data.table(expand.grid(ID2 = ids, ID1 = ids, stringsAsFactors=FALSE))

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
        if(verbose) message("    No clusters identified")
    }

    # add in identity matrix of unrelated samples
    unrel.id <- setdiff(sample.include, names(mem))
    if(verbose) message(length(unrel.id), " samples with no relatives included")

    if(length(unrel.id) > 0) {
        if(is.null(diag.value)){
            # data for the diagonal
            ddat <- x[ID1 == ID2, c('ID1', 'value')][data.table(ID1 = unrel.id), on = 'ID1']
            blocks[[clu$no + 1]] <- Diagonal(n = nrow(ddat), x = ddat$value)
            block.id[[clu$no + 1]] <- ddat$ID1
        }else{
            blocks[[clu$no + 1]] <- Diagonal(n = length(unrel.id), x = diag.value)
            block.id[[clu$no + 1]] <- unrel.id
        }
    }

    # create block diagonal matrix
    if(length(blocks) > 1){
        if(verbose) message("Putting all samples together into one block diagonal matrix")
        mat_sparse <- bdiag(blocks)
    }else{
        mat_sparse <- Matrix(blocks[[1]])
    }

    # ids of samples
    mat.id <- unlist(block.id)
    rownames(mat_sparse) <- mat.id
    colnames(mat_sparse) <- mat.id

    # set the matrix to symmetric to save memory
    mat_sparse <- as(mat_sparse, "symmetricMatrix")
    return(mat_sparse)
}




setGeneric("pcrelateToMatrix", function(pcrelobj, ...) standardGeneric("pcrelateToMatrix"))

setOldClass("pcrelate")
setMethod("pcrelateToMatrix",
          "pcrelate",
          function(pcrelobj, sample.include = NULL, thresh = NULL, scaleKin = 2, verbose = TRUE) {
              .pcrelateToMatrix(pcrelobj = pcrelobj,
                                sample.include = sample.include,
                                thresh = thresh,
                                scaleKin = scaleKin,
                                verbose = verbose)
          })

.pcrelateToMatrix <- function(pcrelobj, sample.include = NULL, thresh = NULL, scaleKin = 2, verbose = TRUE){
    # keep R CMD check from warning about undefined global variables
    f <- ID1 <- ID2 <- kin <- NULL

    # get the diagonals
    x <- as.data.table(pcrelobj$kinSelf)[, c('ID', 'f')]
    setnames(x, 'ID', 'ID1')
    x[, ID2 := ID1]
    x[, kin := 0.5*(1 + f)][, f := NULL]

    # append the off-diagonal
    x <- rbind(x, as.data.table(pcrelobj$kinBtwn)[,c('ID1', 'ID2', 'kin')])

    # scale the values
    x[, kin := scaleKin*kin]

    # set kin name to value
    setnames(x, 'kin', 'value')
    # call makeSparseMatrix
    makeSparseMatrix(x = x, thresh = thresh, sample.include = sample.include, diag.value = NULL, verbose = verbose)
}


setGeneric("kingToMatrix", function(king, ...) standardGeneric("kingToMatrix"))

setMethod("kingToMatrix",
          "character",
          function(king, estimator = c("PropIBD", "Kinship"), sample.include = NULL, thresh = NULL, verbose = TRUE) {
              .kingToMatrix(file.king = king,
                            estimator = estimator,
                            sample.include = sample.include,
                            thresh = thresh,
                            verbose = verbose)
          })

setOldClass("snpgdsIBDClass")
setMethod("kingToMatrix",
          "snpgdsIBDClass",
          function(king, sample.include = NULL, thresh = 2^(-11/2), verbose = TRUE) {
              ID <- king$sample.id
              king <- king$kinship
              dimnames(king) <- list(ID, ID)
              makeSparseMatrix(x = king, thresh = thresh, sample.include = sample.include, diag.value = 0.5, verbose = verbose)
          })

.readKing <- function(x, estimator) {
    cols <- intersect(names(fread(x, nrows=0)), c("ID1", "ID2", estimator))
    if (!(estimator %in% cols)) stop("Column ", estimator, " requested but not present in file")
    fread(x, select=cols, colClasses=list(character=c("ID1", "ID2")))
}

.kingToMatrix <- function(file.king, estimator = c("PropIBD", "Kinship"), sample.include = NULL, thresh = NULL, verbose = TRUE){
    # keep R CMD check from warning about undefined global variables
    PropIBD <- Kinship <- value <- NULL
    #ID1 <- ID2 <- PropIBD <- Kinship <- value <- NULL
    #`.` <- function(...) NULL

    # check argument values
    estimator <- match.arg(estimator)
    if(estimator == 'PropIBD'){
        if(verbose) message('Reading in PropIBD estimates from KING --ibdseg output...')
    }else if(estimator == 'Kinship'){
        if(verbose) message('Reading in Kinship estimates from KING --kinship output...')
    }

    # if multiple input files
    if(length(file.king) > 1){
        king <- lapply(file.king, .readKing, estimator = estimator)
        # pull out column names in common
        cnames <- Reduce(intersect, lapply(king, colnames))
        # subset and rbind
        king <- rbindlist(lapply(king, function(x){ x[, colnames(x) %in% cnames, with = FALSE]} ))

    # one input file
    }else{
        king <- .readKing(file.king, estimator = estimator)
    }

    # subset to needed columns
    if(estimator == 'PropIBD'){
        king <- king[, value := 0.5*PropIBD][, PropIBD := NULL]
    }else if(estimator == 'Kinship'){
        setnames(king, 'Kinship', 'value')
    }

    # check for duplicate pairs
    setkeyv(king, c('ID1', 'ID2'))
    if(any(duplicated(king))){
        stop('Some sample pairs are provided multiple times in file.king; please only provide one value per sample pair')
    }

    # call makeSparseMatrix
    makeSparseMatrix(x = king, thresh = thresh, sample.include = sample.include, diag.value = 0.5, verbose = verbose)
}
