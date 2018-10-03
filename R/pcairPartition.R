pcairPartition <- function(kinobj, divobj,
                           kin.thresh = 2^(-11/2),
                           div.thresh = -2^(-11/2),
                           unrel.set = NULL,
                           sample.include = NULL,
                           verbose = TRUE){

    # sample.id from kinship matrix
    kin.id.all <- .readSampleId(kinobj)
    # sample.id from divergence matrix
    div.id.all <- .readSampleId(divobj)

    # check for ids provided
    if(is.null(kin.id.all) | is.null(div.id.all)) {
        stop('colnames must be provided for kinobj and divobj')
    }

    # logical indicators of which samples are in both sets
    kin.read <- kin.id.all %in% div.id.all
    div.read <- div.id.all %in% kin.id.all
    if(!(all(kin.read) & all(div.read))){
        warning('kinobj and divobj contain non-matching samples; only partitioning those present in both objects')
    }

    # filter to samples in sample.include
    if(!is.null(sample.include)){
        if(!all(sample.include %in% kin.id.all) | !all(sample.include %in% div.id.all)){
            warning('some samples in sample.include are not in kinobj or divobj; they will not be included')
        }

        # subset the read indicators to samples in sample.include
        kin.read <- kin.read & (kin.id.all %in% sample.include)
        div.read <- div.read & (div.id.all %in% sample.include)

        if(verbose) message('Working with ', sum(kin.read), 'samples')
    }

    # checks on unrel.set
    if(!is.null(unrel.set)){
        if(!all(unrel.set %in% kin.id.all) | !all(unrel.set %in% div.id.all)){
            warning('some samples in unrel.set are not in kinobj or divobj; they will not be included')
            # subset unrel.set to only those in kinobj and divobj
            unrel.set <- unrel.set[(unrel.set %in% kin.id.all) & (unrel.set %in% div.id.all)]
        }

        # if using sample.include, only keep unrel.set in sample.include
        if(!is.null(sample.include)){
            unrel.set <- unrel.set[unrel.set %in% sample.include]
        }

        # if no samples left in unrel.set; don't use it downstream
        if(length(unrel.set) == 0){
            unrel.set <- NULL
        }
    }
    

    if(verbose) message('Identifying relatives for each sample using kinship threshold ', kin.thresh)
    # create a vector of ids for samples to be read
    kin.id <- kin.id.all[kin.read]
    # create list of relatives for each sample
    rellist <- .apply(kinobj, MARGIN = 2, 
                      FUN = function(x){ kin.id[x > kin.thresh] },
                      selection = list(kin.read, kin.read))
    names(rellist) <- kin.id
    # remove self from rellist
    for(i in 1:length(rellist)){
        rellist[[i]] <- rellist[[i]][rellist[[i]] != names(rellist)[i]]
    }

    # compute number of relatives for each sample
    nrel <- sapply(rellist, length)
    # set aside samples with no relatives
    unrels <- names(nrel[nrel == 0])
    # subset rellist, nrel to those with relatives
    rellist <- rellist[!(names(rellist) %in% unrels)]
    nrel <- nrel[nrel > 0]


    # update logicial indicator of which samples in kinobj still need to be read
    kin.read <- kin.id.all %in% names(nrel)
    # update the vector of ids we are reading
    kin.id <- kin.id.all[kin.read]
    # compute 'total kinship' values
    kinsum <- .apply(kinobj, MARGIN = 2, 
                     FUN = function(x){ sum(x[x > kin.thresh]) },
                     selection = list(kin.read, kin.read))
    names(kinsum) <- kin.id
    kinsum <- unlist(kinsum)


    if(verbose) message('Identifying pairs of divergent samples using divergence threshold ', div.thresh)
    # logical indicator of which samples in divobj need divergence measures
    div.read.col <- div.id.all %in% kin.id
    # vector of ids we are reading
    div.id <- div.id.all[div.read]
    div.id.col <- div.id.all[div.read.col]
    # create list of divergent pairs for each sample
    divlist <- .apply(divobj, MARGIN = 2, 
                      FUN = function(x){ div.id[x < div.thresh] }, 
                      selection = list(div.read, div.read.col))
    names(divlist) <- div.id.col

    # create a vector matching ids of rellist and divlist
    idx <- match(names(divlist), names(rellist))
    # not divergent if actually related
    for(i in 1:length(divlist)){
        j <- idx[i]
        divlist[[i]] <- divlist[[i]][!(divlist[[i]] %in% rellist[[j]])]
    }

    # compute number of divergent pairs for each sample
    ndiv <- sapply(divlist, length)

    # clean up
    rm(divlist); rm(div.id.all); rm(div.read); rm(div.read.col); rm(div.id); rm(div.id.col)


    # empty vector to store related set
    rels <- NULL


    # take care of user specified unrel.set
    if(!is.null(unrel.set)){
        if(verbose) message('Forcing samles specified in unrel.set into the unrelated set')
        # add these samples to unrels
        unrels <- unique(append(unrels, unrel.set))

        # identify samples related to a sample in unrel.set
        rel.new <- names(rellist)[sapply(rellist, function(x){ any(x %in% unrel.set) })]
        # append rel.new to the master relative list
        rels <- append(rels, rel.new)

        # remove unrel.set and rel.new from rellist
        keep <- !(names(nrel) %in% c(unrel.set, rel.new))
        rellist <- rellist[keep]
        # recompute nrel
        nrel <- sapply(rellist, length)

        # identify any with no relatives left
        unrel.new <- names(nrel[nrel == 0])
        if(length(unrel.new) > 0){
            # append unrel.new to the master unrelated list
            unrels <- append(unrels, unrel.new)
            
            # remove unrel.new from rellist and nrel
            keep <- !(names(nrel) %in% unrel.new)
            rellist <- rellist[keep]
            nrel <- nrel[keep]
        }
    }


    if(verbose) message('Paritioning samples into unrelated and related sets...')
    # iterate
    iter <- 0
    while(length(nrel) > 0){
        iter <- iter + 1
        if(verbose & iter %% 1000 == 0) message('    ...', iter, ' samples added to related.set...')
        
        # who has the most relatives
        rel.new <- names(nrel)[nrel == max(nrel)]
        if(length(rel.new) > 1){
            # who has the least divergent pairs
            ndiv2 <- ndiv[names(ndiv) %in% rel.new]
            rel.new <- names(ndiv2)[ndiv2 == min(ndiv2)]
            if(length(rel.new) > 1){
                # who has the lowest 'total kinship'
                kinsum2 <- kinsum[names(kinsum) %in% rel.new]
                rel.new <- names(kinsum2)[kinsum2 == min(kinsum2)][1]
            }
        }

        # append rel.new to the master relative list
        rels <- append(rels, rel.new)
        
        # remove rel.new from rellist and nrel
        keep <- names(nrel) != rel.new
        rellist <- rellist[keep]
        nrel <- nrel[keep]

        # subtract 1 from nrel for those related to rel.new
        nrel <- nrel - ifelse(sapply(rellist, function(x){ rel.new %in% x }), 1, 0)
        
        # identify any with no relatives left
        unrel.new <- names(nrel[nrel == 0])
        if(length(unrel.new) > 0){
            # append unrel.new to the master unrelated list
            unrels <- append(unrels, unrel.new)
            
            # remove unrel.new from rellist and nrel
            keep <- !(names(nrel) %in% unrel.new)
            rellist <- rellist[keep]
            nrel <- nrel[keep]
        }
    }

    # return results
    return(list(rels = rels, unrels = unrels))

}


.readSampleId <- function(x, ...) UseMethod(".readSampleId", x)
.readSampleId.matrix <- function(x) {
    colnames(x)
}

.readSampleId.Matrix <- function(x) {
    colnames(x)
}

.readSampleId.gds.class <- function(x) {
    read.gdsn(index.gdsn(x, 'sample.id'))
}


.apply <- function(x, MARGIN, FUN, selection, ...) UseMethod(".apply", x)
.apply.matrix <- function(x, MARGIN, FUN, selection) {
    apply(x[selection[[1]], selection[[2]]], MARGIN = MARGIN, FUN = FUN)
}

.apply.Matrix <- function(x, MARGIN, FUN, selection) {
    x <- x[selection[[1]], selection[[2]]]
    ans <- list()
    if (MARGIN == 1) {
        for (i in 1:nrow(x)) {
            ans[[i]] <- FUN(x[i,])
        }
        names(ans) <- rownames(x)
    } else if (MARGIN == 2) {
        for (i in 1:ncol(x)) {
            ans[[i]] <- FUN(x[,i])
        }
        names(ans) <- colnames(x)
    } else {
        stop("MARGIN must be 1 or 2")
    }
    simplify2array(ans)
}

.apply.gds.class <- function(x, MARGIN, FUN, selection) {
    apply.gdsn(index.gdsn(x, 'kinship'), margin = MARGIN, FUN = FUN,
               selection = selection)
}


.apply.Matrix_v2 <- function (x, MARGIN, FUN, selection, maxelem = 2^31 - 1){
    
    # subset to selection
    x <- x[selection[[1]], selection[[2]]]
    
    # determine number of blocks needed
    nr <- as.numeric(nrow(x))
    nc <- as.numeric(ncol(x))
    nblock <- ceiling(nr*nc/maxelem)

    if(nblock > 1){
        
        if(MARGIN  == 1){
            # create blocks
            bsize <- floor(nr/nblock)
            if(bsize < 2){
                stop('block size is too small; increase maxelem')
            }
            start <- 0:(nblock-1)*bsize + 1
            stop <- c(1:(nblock-1)*bsize, nr)
            
            # loop through blocks
            ans <- list()
            for(i in 1:nblock){
                ans[[i]] <- apply(x[start[i]:stop[i], ], 1, FUN)
            }
            
        }else if(MARGIN == 2){
            # create blocks
            bsize <- floor(nc/nblock)   
            start <- 0:(nblock-1)*bsize + 1
            stop <- c(1:(nblock-1)*bsize, nc)
            
            # loop through blocks
            ans <- list()
            for(i in 1:nblock){
                ans[[i]] <- apply(x[ ,start[i]:stop[i]], 2, FUN)
            }
            
        }else {
            stop("MARGIN must be 1 or 2")
        }
            
        # unlist the top level
        if(length(unique(sapply(ans, class))) > 1){
            stop('not all elements in the output of .apply.Matrix have the same class')
        }
        unlist(ans, recursive = FALSE)  
            
    }else{
        apply(x, MARGIN, FUN)
    }
}
