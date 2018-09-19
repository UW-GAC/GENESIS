kingToMatrix <- function(file.king, sample.include = NULL, thresh = NULL, verbose = TRUE){
    
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
        king <- king[,.(ID1, ID2, PropIBD)]
        king <- king[, Kinship := 0.5*PropIBD][, PropIBD := NULL]
    }else if('Kinship' %in% colnames(king)){
        if(verbose) message('Inferred to be KING --robust output')
        king <- king[,.(ID1, ID2, Kinship)]
    }else{
        stop('All files in file.king must contain a column called `PropIBD` (KING --ibdseg output) or `Kinship` (KING --robust output)')
    }
        
    # check for duplicate pairs
    setkey(king, ID1, ID2)
    if(any(duplicated(king))){
        stop('Some sample pairs are provided multiple times in file.king; please only provide one value per sample pair')
    }
        
    if(!is.null(sample.include)){
        # subset king data to samples in sample.include
        if(verbose) message('Using ', length(sample.include), ' samples in sample.include')
        king <- king[ID1 %in% sample.include & ID2 %in% sample.include]
    }else{
        # get list of all samples in the file
        sample.include <- sort(unique(c(king$ID1, king$ID2)))
        if(verbose) message('Using ', length(sample.include), ' samples in file.king')
    }
    
    if(verbose) message('Identifying clusters of relatives...')
    # create linked graph (each of the first two columns specifies edges)
    if(is.null(thresh)){
        g <- igraph::graph_from_data_frame(king)
    }else{
        g <- igraph::graph_from_data_frame(king[Kinship > thresh])
    }
    # extract cluster membership
    clu <- igraph::components(g)
    mem <- clu$membership
    if(verbose) message('    ', length(mem), ' relatives in ', clu$no, ' clusters')
    
    
    # make an identity matrix of unrelated samples
    unrel.id <- setdiff(sample.include, names(mem))
    if(verbose) message(length(unrel.id), ' samples with no relatives included...')
    mat_unrels <- 0.5*Diagonal(length(unrel.id))
    
    
    if(clu$no > 0){
        if(verbose) message('Creating block matrices for clusters...')
        blocks <- list()
        block.id <- list()
        for(i in 1:clu$no){
            # samples in the cluster
            ids <- names(mem[mem == i])
            # create a table for all pairs in the cluster
            allpairs <- as.data.table(expand.grid(ID2 = ids, ID1 = ids))

            # merge
            sub <- king[ID1 %in% ids & ID2 %in% ids][allpairs, on = c('ID1', 'ID2')]
            # set pairs not included in KING data to 0
            sub[is.na(Kinship), Kinship := 0]

            # cast to a matrix
            submat <- reshape2::acast(data = sub, formula = ID1 ~ ID2, value.var = 'Kinship')
            # put the values on both sides of the diagonal
            submat <- submat + t(submat)
            # set the diagonal to 0.5
            diag(submat) <- 0.5

            # store in the list
            blocks[[i]] <- submat
            block.id[[i]] <- rownames(submat)
        }
        # turn the list into a block-diagonal matrix
        mat_rels <- bdiag(blocks)
        # ids of samples with relatives
        rel.id <- unlist(block.id)
        
        if(verbose) message('Putting all samples together into one block diagonal matrix')
        # create a block diagonal matrix with all samples
        mat_sparse <- bdiag(mat_rels, mat_unrels)
        rownames(mat_sparse) <- c(rel.id, unrel.id)
        colnames(mat_sparse) <- c(rel.id, unrel.id)
        
        return(mat_sparse)
        
    }else{
        rownames(mat_unrels) <- unrel.id
        colnames(mat_unrels) <- unrel.id
        return(mat_unrels)
    }
}
