setGeneric("pcrelate", function(gdsobj, ...) standardGeneric("pcrelate"))

setMethod("pcrelate",
          "GenotypeData",
          function(gdsobj, ...) {
              .Defunct("pcrelate", msg="The GenotypeData method for pcrelate is defunct. Use the GenotypeIterator method instead.")
          })

setMethod("pcrelate",
          "SeqVarData",
          function(gdsobj, ...) {
              .Defunct("pcrelate", msg="The SeqVarData method for pcrelate is defunct. Use the SeqVarIterator method instead.")
          })

setMethod("pcrelate",
          "GenotypeIterator",
          function(gdsobj,
                   pcs,
                   scale = c('overall', 'variant', 'none'),
                   ibd.probs = TRUE,
                   sample.include = NULL,
                   training.set = NULL,
                   sample.block.size = 5000,
                   maf.thresh = 0.01,
                   maf.bound.method = c('filter', 'truncate'),
                   small.samp.correct = TRUE,
                   #num.cores = 1,
                   verbose = TRUE) {
              .pcrelate(gdsobj, 
                        pcs = pcs,
                        scale = scale,
                        ibd.probs = ibd.probs,
                        sample.include = sample.include,
                        training.set = training.set,
                        sample.block.size = sample.block.size,
                        maf.thresh = maf.thresh,
                        maf.bound.method = maf.bound.method,
                        small.samp.correct = small.samp.correct,
                        #num.cores = num.cores,
                        verbose = verbose)
          })

setMethod("pcrelate",
          "SeqVarIterator",
          function(gdsobj,
                   pcs,
                   scale = c('overall', 'variant', 'none'),
                   ibd.probs = TRUE,
                   sample.include = NULL,
                   training.set = NULL,
                   sample.block.size = 5000,
                   maf.thresh = 0.01,
                   maf.bound.method = c('filter', 'truncate'),
                   small.samp.correct = TRUE,
                   #num.cores = 1,
                   verbose = TRUE) {
              filt <- seqGetFilter(gdsobj)
              out <- .pcrelate(gdsobj, 
                               pcs = pcs,
                               scale = scale,
                               ibd.probs = ibd.probs,
                               sample.include = sample.include,
                               training.set = training.set,
                               sample.block.size = sample.block.size,
                               maf.thresh = maf.thresh,
                               maf.bound.method = maf.bound.method,
                               small.samp.correct = small.samp.correct,
                               #num.cores = num.cores,
                               verbose = verbose)
              seqSetFilter(gdsobj,
                           sample.sel=filt$sample.sel,
                           variant.sel=filt$variant.sel,
                           verbose=FALSE)
              out
          })

.pcrelate <- function(gdsobj,
                      pcs,
                      scale = c('overall', 'variant', 'none'),
                      ibd.probs = TRUE,
                      sample.include = NULL,
                      training.set = NULL,
                      sample.block.size = 5000,
                      maf.thresh = 0.01,
                      maf.bound.method = c('filter', 'truncate'),
                      small.samp.correct = TRUE,
                      #num.cores = 1,
                      verbose = TRUE){

    # checks
    scale <- match.arg(scale)
    maf.bound.method <- match.arg(maf.bound.method)
    sample.include <- samplesGdsOrder(gdsobj, sample.include)
    .pcrelateChecks(pcs = pcs, scale = scale, ibd.probs = ibd.probs, sample.include = sample.include, training.set = training.set, 
                    maf.thresh = maf.thresh)
    
    # set up number of cores
    ## sys.cores <- parallel::detectCores(logical = TRUE)
    ## doMC::registerDoMC(cores = min(c(num.cores, sys.cores)))
    ## if(verbose) message('Using ', min(c(num.cores, sys.cores)), ' CPU cores')

    # number of sample blocks
    nsampblock <- ceiling(length(sample.include)/sample.block.size)

    # check for small sample correction
    if(small.samp.correct){
        small.samp.correct <- (nsampblock == 1) & (scale != 'none')
        if(!small.samp.correct) warning('small.samp.correct can only be used when all samples are analyzed in one block and `scale != none`')
    }

    # list of samples in each block
    if(nsampblock > 1){
        samp.blocks <- unname(split(sample.include, cut(1:length(sample.include), nsampblock)))
    }else{
        samp.blocks <- list(sample.include)
    }

    if(nsampblock == 1){
        if(verbose) message(length(sample.include), ' samples to be included in the analysis...')

        # create matrix of PCs
        V <- .createPCMatrix(pcs = pcs, sample.include = sample.include)

        # matrix product of V
        VVtVi <- .calcISAFBetaPCProd(V = V, training.set = training.set, verbose = verbose)

        snp.blocks <- .snpBlocks(gdsobj)
        nsnpblock <- length(snp.blocks)
        if(verbose) message('Running PC-Relate analysis for ', length(sample.include), ' samples using ', length(unlist(snp.blocks)), ' SNPs in ', nsnpblock, ' blocks...')

        # for each snp block
        matList <- foreach(k = 1:nsnpblock, .combine = .matListCombine, .inorder = FALSE, .multicombine = FALSE) %do% {
            if(verbose) message('    Running block ', k, '...')
            # read genotype data for the block
            G <- .readGeno(gdsobj, sample.include, snp.index = snp.blocks[[k]])

            # calculate ISAF betas
            beta <- .calcISAFBeta(G = G, VVtVi = VVtVi)

            # calculate PC-Relate estimates
            .pcrelateVarBlock(G = G, beta = beta, V = V, idx = 1:nrow(V), jdx = 1:nrow(V), scale = scale, ibd.probs = ibd.probs, maf.thresh = maf.thresh, maf.bound.method = maf.bound.method)
        }

        # take ratios to compute final estimates
        estList <- .pcrelateCalcRatio(matList = matList, scale = scale, ibd.probs = ibd.probs)
        rm(matList)

        # cast to data.tables
        kinSelf <- data.table(ID = rownames(estList$kin), f = 2*diag(estList$kin) - 1, nsnp = diag(estList$nsnp))
        kinBtwn <- .estListToDT(estList, drop.lower = TRUE)
        rm(estList)


    }else if(nsampblock > 1){
        if(verbose) message(length(sample.include), ' samples to be included in the analysis, split into ', nsampblock, ' blocks...')

        # calculate betas for individual specific allele frequencies
        beta <- calcISAFBeta(gdsobj = gdsobj, pcs = pcs, sample.include = sample.include, training.set = training.set, verbose = verbose)		
        ### beta is a matrix of variants x PCs; needs to be saved, not sure of the best format or the best way to split this up ###
        
        # compute estimates for current (pair of) sample block(s)
        ### this is where we would parallelize with multiple jobs ###
        kinSelf <- NULL
        kinBtwn <- NULL
        for(i in 1:nsampblock){
            for(j in i:nsampblock){
                if(verbose) message('Running PC-Relate analysis for sample block pair (', i, ',', j, ')')
                # compute estimates for the (pair of) sample block(s)
                tmp <- pcrelateSampBlock(	gdsobj = gdsobj, pcs = pcs, betaobj = beta, sample.include.block1 = samp.blocks[[i]], sample.include.block2 = samp.blocks[[j]],
                                         scale = scale, ibd.probs = ibd.probs,
                                         maf.thresh = maf.thresh, maf.bound.method = maf.bound.method, verbose = verbose)

                # update results with this (pair of) sample block(s)
                if(i == j) kinSelf <- rbind(kinSelf, tmp$kinSelf)
                kinBtwn <- rbind(kinBtwn, tmp$kinBtwn)
            }
        }
    }

    ### post processing after putting all samples together

    # order samples
    setkeyv(kinSelf, 'ID')
    setkeyv(kinBtwn, c('ID1', 'ID2'))

    # correct kinship - small sample
    if(small.samp.correct){
        if(verbose) message('Performing Small Sample Correction...')
        out <- correctKin(kinBtwn = kinBtwn, kinSelf = kinSelf, pcs = pcs, sample.include = sample.include)
        kinBtwn <- out$kinBtwn
        kinSelf <- out$kinSelf
    }

    # correct k2 - HW departure and small sample
    if(ibd.probs) kinBtwn <- correctK2(kinBtwn = kinBtwn, kinSelf = kinSelf, small.samp.correct = small.samp.correct, pcs = pcs, sample.include = sample.include)

    # use alternate k0 estimator for non-1st degree relatives
    if(ibd.probs) kinBtwn <- correctK0(kinBtwn = kinBtwn)
    
    # return output
    out <- list(kinBtwn = as.data.frame(kinBtwn), kinSelf = as.data.frame(kinSelf))
    class(out) <- "pcrelate"
    return(out)
}




.pcrelateChecks <- function(pcs, scale, ibd.probs, sample.include, training.set, maf.thresh){
    # check parameters
    if(scale == 'none' & ibd.probs) stop('`ibd.probs` must be FALSE when `scale` == none')
    if(maf.thresh < 0 | maf.thresh > 0.5) stop('maf.thresh must be in [0,0.5]')
    # check training.set
    if(!is.null(training.set) & !all(training.set %in% sample.include)) stop('All samples in training.set must be in sample.include')
    # check pcs
    if(!is.matrix(pcs) | is.null(rownames(pcs))) stop('pcs should be a matrix of PCs with rownames set to sample.ids')
    if(!all(sample.include %in% rownames(pcs))) stop('All samples in sample.include must be in pcs')
}


### get sample ids in same order as gdsobj
samplesGdsOrder <- function(gdsobj, sample.include) {
    sample.id <- .readSampleId(gdsobj)
    if (!is.null(sample.include)) {
        sample.id <- intersect(sample.id, sample.include)
    }
    return(as.character(sample.id))
}


### function to match samples and create PC matrix
.createPCMatrix <- function(pcs, sample.include){
    # subset and re-order pcs if needed
    V <- pcs[match(sample.include, rownames(pcs)), , drop = FALSE]
    # append intercept
    V <- cbind(rep(1, nrow(V)), V)
    return(V)
}


# function to get 
.calcISAFBetaPCProd <- function(V, training.set, verbose = TRUE){
    if(!is.null(training.set)){
        idx <- rownames(V) %in% training.set
        VVtVi <- tcrossprod(V[idx,], chol2inv(chol(crossprod(V[idx,]))))
        if(verbose) message('Betas for ', ncol(V) - 1, ' PC(s) will be calculated using ', sum(idx), ' samples in training.set...')
    }else{
        idx <- NULL
        VVtVi <- tcrossprod(V, chol2inv(chol(crossprod(V))))
        if(verbose) message('Betas for ', ncol(V) - 1, ' PC(s) will be calculated using all ', nrow(V), ' samples in sample.include...')
    }
    return(list(val = VVtVi, idx = idx))
}

# function to do actual calculation of betas
.calcISAFBeta <- function(G, VVtVi){
    # impute missing genotype values
    G <- .meanImpute(G, freq=0.5*colMeans(G, na.rm=TRUE))

    # calculate beta
    if(is.null(VVtVi$idx)){
        if(!identical(rownames(G), rownames(VVtVi$val))) stop('sample order in genotypes and pcs do not match')
        beta <- crossprod(G, VVtVi$val)
    }else{
        if(!identical(rownames(G)[VVtVi$idx], rownames(VVtVi$val))) stop('sample order in genotypes and pcs do not match')
        beta <- crossprod(G[VVtVi$idx, ], VVtVi$val)
    }
    return(beta)
}



### this function does the pcrelate estimation for a variant block
.pcrelateVarBlock <- function(G, beta, V, idx, jdx, scale, ibd.probs, maf.thresh, maf.bound.method){

    # make sure order of G, beta, and V all line up
    if(!identical(colnames(G), rownames(beta))) stop('G and beta do not match')
    if(!identical(rownames(G), rownames(V))) stop('G and V do not match')

    # estimate individual specific allele frequencies
    mu <- .estISAF(beta = beta, V = V, bound.thresh = maf.thresh, bound.method = maf.bound.method)

    # compute values for estimates
    if(scale == 'overall'){
        matList <- .pcrCalcOvr(G = G, mu = mu, idx = idx, jdx = jdx, ibd.probs = ibd.probs)
    }else if(scale == 'variant'){
        matList <- .pcrCalcVar(G = G, mu = mu, idx = idx, jdx = jdx, ibd.probs = ibd.probs)
    }else if(scale == 'none'){
        matList <- .pcrCalcNone(G = G, mu = mu, idx = idx, jdx = jdx)
    }

    return(matList)
}

.matListCombine <- function(...){
    mapply(FUN = "+", ..., SIMPLIFY = FALSE)
}


### these functions are for estimating individual specific allele frequencies
.estISAF <- function(beta, V, bound.thresh, bound.method){
    # get ISAF estimates (i.e. 0.5*fitted values)
    mu <- 0.5*tcrossprod(V, beta)
    # fix boundary cases
    if(bound.method == 'truncate'){
        mu <- apply(mu, 2, .muBoundaryTrunc, thresh = bound.thresh)
    }else if(bound.method == 'filter'){
        mu <- apply(mu, 2, .muBoundaryFilt, thresh = bound.thresh)
    }else{
        stop('bound.method must be one of `truncate` or `filter`')
    }
    return(mu)
}

# function to truncate boundary issues with ISAFs
.muBoundaryTrunc <- function(x, thresh){
    x[x < thresh] <- thresh
    x[x > (1 - thresh)] <- (1 - thresh)
    return(x)
}

# function to filter boundary issues with ISAFs
.muBoundaryFilt <- function(x, thresh){
    x[x < thresh] <- NA
    x[x > (1 - thresh)] <- NA
    return(x)
}



### all of the functions below are used for computing kinship and ibd estimates ###
.pcrCalcOvr <- function(G, mu, ibd.probs, idx, jdx){
    # compute mu(1 - mu)
    muqu <- mu*(1 - mu)

    # compute number of observed snps by pair
    nonmiss <- !(is.na(G) | is.na(mu))
    nsnp <- tcrossprod(nonmiss[idx,,drop=F], nonmiss[jdx,,drop=F])

    # index of missing values
    filt.idx <- which(!nonmiss)

    # compute kinship values
    kinList <- .pcrCalcKinOvr(G, mu, muqu, filt.idx, idx, jdx)

    if(ibd.probs){
        # compute ibd values
        ibdList <- .pcrCalcIBDOvr(G, mu, muqu, nonmiss, filt.idx, idx, jdx)

        return(list(kinNum = kinList$kinNum, 
                    kinDen = kinList$kinDen,
                    k0Num = ibdList$k0Num,
                    k0Den = ibdList$k0Den,
                    k2Num = ibdList$k2Num,
                    k2Den = ibdList$k2Den,
                    nsnp = nsnp))

    }else{
        return(list(kinNum = kinList$kinNum,
                    kinDen = kinList$kinDen,
                    nsnp = nsnp))
    }
}

.pcrCalcKinOvr <- function(G, mu, muqu, filt.idx, idx, jdx){
    # residuals
    R <- G - 2*mu

    # set missing values to 0 (i.e. no contribution from that variant)
    R[filt.idx] <- 0
    muqu[filt.idx] <- 0

    # numerator (crossprod of residuals)
    kinNum <- tcrossprod(R[idx,,drop=F], R[jdx,,drop=F])
    # denominator
    kinDen <- tcrossprod(sqrt(muqu[idx,,drop=F]), sqrt(muqu[jdx,,drop=F]))

    return(list(kinNum = kinNum, kinDen = kinDen))
}

.pcrCalcIBDOvr <- function(G, mu, muqu, nonmiss, filt.idx, idx, jdx){
    # opposite allele frequency
    qu <- 1 - mu

    # indicator matrices of homozygotes
    Iaa <- G == 0 & nonmiss
    IAA <- G == 2 & nonmiss

    # dominance coded genotype matrix
    Gd <- mu
    Gd[G == 1 & nonmiss] <- 0
    Gd[IAA] <- qu[IAA]
    Gd <- Gd - muqu

    # set missing values to 0 (i.e. no contribution from that variant)
    Iaa[filt.idx] <- 0
    IAA[filt.idx] <- 0
    Gd[filt.idx] <- 0
    mu[filt.idx] <- 0
    qu[filt.idx] <- 0
    muqu[filt.idx] <- 0

    # k0	
    k0Num <- tcrossprod(IAA[idx,,drop=F], Iaa[jdx,,drop=F]) + tcrossprod(Iaa[idx,,drop=F], IAA[jdx,,drop=F])
    k0Den <- tcrossprod(mu[idx,,drop=F]^2, qu[jdx,,drop=F]^2) + tcrossprod(qu[idx,,drop=F]^2, mu[jdx,,drop=F]^2)

    # k2
    k2Num <- tcrossprod(Gd[idx,,drop=F], Gd[jdx,,drop=F])
    k2Den <- tcrossprod(muqu[idx,,drop=F], muqu[jdx,,drop=F])

    return(list(k0Num = k0Num, k0Den = k0Den, k2Num = k2Num, k2Den = k2Den))
}


.pcrCalcVar <- function(G, mu, ibd.probs, idx, jdx){
    # compute mu(1 - mu)
    muqu <- mu*(1 - mu)

    # compute number of observed snps by pair
    nonmiss <- !(is.na(G) | is.na(mu))
    nsnp <- tcrossprod(nonmiss[idx,,drop=F], nonmiss[jdx,,drop=F])

    # index of missing values
    filt.idx <- which(!nonmiss)

    # compute kinship values
    kinNum <- .pcrCalcKinVar(G, mu, muqu, filt.idx, idx, jdx)

    if(ibd.probs){
        # compute ibd values
        ibdList <- .pcrCalcIBDVar(G, mu, muqu, nonmiss, filt.idx, idx, jdx)

        return(list(kinNum = kinNum,
                    k0Num = ibdList$k0Num,
                    k2Num = ibdList$k2Num,
                    nsnp = nsnp))

    }else{
        return(list(kinNum = kinNum, nsnp = nsnp))
    }
}

.pcrCalcKinVar <- function(G, mu, muqu, filt.idx, idx, jdx){
    # scaled residuals
    R <- (G - 2*mu)/sqrt(muqu)

    # set missing values to 0 (i.e. no contribution from that variant)
    R[filt.idx] <- 0

    # numerator (crossprod of scaled residuals)
    kinNum <- tcrossprod(R[idx,,drop=F], R[jdx,,drop=F])
    
    return(kinNum)
}

.pcrCalcIBDVar <- function(G, mu, muqu, nonmiss, filt.idx, idx, jdx){
    # opposite allele frequency
    qu <- 1 - mu

    # indicator matrices of homozygotes
    Iaa <- G == 0 & nonmiss
    IAA <- G == 2 & nonmiss

    # dominance coded genotype matrix
    Gd <- mu
    Gd[G == 1 & nonmiss] <- 0
    Gd[IAA] <- qu[IAA] 

    # scale 
    Iaa <- 0.5*Iaa/qu^2 # factor of 1/2 for IAA*Iaa comes in from splitting up the two terms AA,aa vs aa,AA
    IAA <- IAA/mu^2
    Gd <- Gd/muqu - 1

    # set missing values to 0 (i.e. no contribution from that variant)
    Iaa[filt.idx] <- 0
    IAA[filt.idx] <- 0
    Gd[filt.idx] <- 0

    # k0
    k0Num <- tcrossprod(IAA[idx,,drop=F], Iaa[jdx,,drop=F]) + tcrossprod(Iaa[idx,,drop=F], IAA[jdx,,drop=F])

    # k2
    k2Num <- tcrossprod(Gd[idx,,drop=F], Gd[jdx,,drop=F])

    return(list(k0Num = k0Num, k2Num = k2Num))
}


.pcrCalcNone <- function(G, mu, idx, jdx){
    # compute number of observed snps by pair
    nonmiss <- !(is.na(G) | is.na(mu))
    nsnp <- tcrossprod(nonmiss[idx,,drop=F], nonmiss[jdx,,drop=F])

    # index of missing values
    filt.idx <- which(!nonmiss)
    rm(nonmiss)

    # compute kinship values
    kin <- .pcrCalcKinNone(G, mu, filt.idx, idx, jdx)
    return(list(kin = kin, nsnp = nsnp))
}

.pcrCalcKinNone <- function(G, mu, filt.idx, idx, jdx){
    # residuals
    R <- G - 2*mu

    # set missing values to 0 (i.e. no contribution from that variant)
    R[filt.idx] <- 0

    # crossprod of scaled residuals
    kin <- tcrossprod(R[idx,,drop=F], R[jdx,,drop=F])

    return(kin)
}


### functions for post processing on variant block level
.pcrelateCalcRatio <- function(matList, scale, ibd.probs){
    # compute final estimates
    if(scale == 'overall'){
        kin <- matList$kinNum/(4*matList$kinDen)
        if(ibd.probs){
            k2 <- matList$k2Num/matList$k2Den
            k0 <- matList$k0Num/matList$k0Den
            return(list(kin = kin, k0 = k0, k2 = k2, nsnp = matList$nsnp))
        }else{
            return(list(kin = kin, nsnp = matList$nsnp))
        }

    }else{
        kin <- matList$kin/(4*matList$nsnp)
        if(ibd.probs){
            k2 <- matList$k2/matList$nsnp
            k0 <- matList$k0/matList$nsnp
            return(list(kin = kin, k0 = k0, k2 = k2, nsnp = matList$nsnp))
        }else{
            return(list(kin = kin, nsnp = matList$nsnp))
        }
    }
}


### functions for post processing on sample block level
.estListToDT <- function(x, drop.lower){
    estDT <- lapply(x, meltMatrix, drop.lower = drop.lower, drop.diag = drop.lower)
    for(k in 1:length(estDT)){
        setnames(estDT[[k]], 'value', names(estDT)[k])
    }
    # merge those data.tables into one data.table
    estDT <- Reduce(merge, estDT)
    return(estDT)
}


### functions for final processing
correctKin <- function(kinBtwn, kinSelf, pcs, sample.include = NULL){
    # keep R CMD check from warning about undefined global variables
    f <- ID1 <- ID2 <- kin <- newval <- value <- NULL
    
    # temporary data.table to store values
    tmp <- kinSelf[, c('ID', 'f')]
    setnames(tmp, c('ID','f'), c('ID1', 'kin'))
    tmp[, ID2 := ID1]
    tmp <- rbind(kinBtwn[, c('ID1', 'ID2', 'kin')], tmp)
    setnames(tmp, 'kin', 'newval')
    setkeyv(tmp, c('ID1', 'ID2'))
    
    # get the PC matrix
    V <- .createPCMatrix(pcs = pcs, sample.include = sample.include)
    
    # make adjustment for each PC
    for(k in 2:ncol(V)){
        Acov <- tcrossprod(V[,k])
        rownames(Acov) <- rownames(V)
        colnames(Acov) <- rownames(V)
        Avec <- meltMatrix(Acov, drop.lower = TRUE, drop.diag = FALSE)
        tmp <- Avec[tmp, on = c('ID1', 'ID2')]
        coef <- lm(formula = as.formula(newval ~ value), data = tmp[newval < 2^(-11/2)])$coef
        tmp[, newval := newval - coef[1] - coef[2]*value]
        tmp[, value := NULL]
    }
    
    # merge back into kinBtwn
    kinBtwn <- tmp[kinBtwn, on = c('ID1', 'ID2')]
    kinBtwn[, kin := newval][, newval := NULL]
    
    # merge back into kinSelf
    tmp <- tmp[ID1 == ID2][, ID2 := NULL]
    setnames(tmp, 'ID1', 'ID')
    kinSelf <- tmp[kinSelf, on = 'ID']
    kinSelf[, f := newval][, newval := NULL]

    return(list(kinBtwn = kinBtwn, kinSelf = kinSelf))
}

correctK2 <- function(kinBtwn, kinSelf, pcs, sample.include = NULL, small.samp.correct = TRUE){
    # keep R CMD check from warning about undefined global variables
    f.1 <- f.2 <- kin <- k2 <- newval <- value <- NULL
    
    # correct k2 for HW departure
    kinBtwn <- merge(kinBtwn, kinSelf[, c('ID', 'f')], by.x = 'ID2', by.y = 'ID')
    setnames(kinBtwn, 'f', 'f.2')
    kinBtwn <- merge(kinBtwn, kinSelf[, c('ID', 'f')], by.x = 'ID1', by.y = 'ID')
    setnames(kinBtwn, 'f', 'f.1')
    kinBtwn[, k2 := k2 - f.1*f.2][, `:=`(f.1 = NULL, f.2 = NULL)]

    if(small.samp.correct){
        # temporary data.table to store values
        tmp <- kinBtwn[, c('ID1', 'ID2', 'kin', 'k2')]
        setnames(tmp, 'k2', 'newval')
        setkeyv(tmp, c('ID1', 'ID2'))

        # get the PC matrix
        V <- .createPCMatrix(pcs = pcs, sample.include = sample.include)

        # make adjustment for each PC
        for(k in 2:ncol(V)){
            Acov <- tcrossprod(V[,k])
            rownames(Acov) <- rownames(V)
            colnames(Acov) <- rownames(V)
            Avec <- meltMatrix(Acov, drop.lower = TRUE, drop.diag = FALSE)
            tmp <- Avec[tmp, on = c('ID1', 'ID2')]
            coef <- lm(formula = as.formula(newval ~ value + I(value^2)), data = tmp[kin < 2^(-11/2)])$coef
            tmp[, newval := newval - coef[1] - coef[2]*value - coef[3]*value^2]
            tmp[, value := NULL]
        }
        
        # make adjustment for kinship
        coef <- lm(formula = as.formula(newval ~ kin), data = tmp[newval < 2^(-9/2)])$coef
        tmp[, newval := newval - coef[1] - coef[2]*kin]
        tmp[, kin := NULL]
        
        # merge back into kinBtwn
        kinBtwn <- tmp[kinBtwn, on = c('ID1', 'ID2')]
        kinBtwn[!is.na(newval), k2 := newval][, newval := NULL]
    }
    
    return(kinBtwn)
}

correctK0 <- function(kinBtwn){
    # keep R CMD check from warning about undefined global variables
    kin <- k0 <- k2 <- NULL
    
    # use alternate k0 estimator for non-1st degree relatives
    kinBtwn[kin < 2^(-5/2), k0 := 1 - 4*kin + k2]
    
    return(kinBtwn)
}


setGeneric("meltMatrix", function(x, ...) standardGeneric("meltMatrix"))

setMethod("meltMatrix",
          "matrix",
          function(x, drop.lower = FALSE, drop.diag = FALSE){
            ID1 <- ID2 <- NULL
            if(drop.lower){
                x[lower.tri(x, diag = drop.diag)] <- NA
            }
            x <- as.data.table(reshape2::melt(x, varnames = c('ID1', 'ID2'), na.rm = TRUE, as.is = TRUE))
            x <- x[,`:=`(ID1 = as.character(ID1), ID2 = as.character(ID2))]
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





### exported function for computing PC betas for individual specific allele frequency calculations ###
calcISAFBeta <- function(gdsobj, pcs, sample.include, training.set = NULL, verbose = TRUE){
    # checks - add some

    # create matrix of PCs
    V <- .createPCMatrix(pcs = pcs, sample.include = sample.include)

    # matrix product of V
    VVtVi <- .calcISAFBetaPCProd(V = V, training.set = training.set, verbose = verbose)
    
    snp.blocks <- .snpBlocks(gdsobj)
    nsnpblock <- length(snp.blocks)
    if(verbose) message('Calculating Indivdiual-Specific Allele Frequency betas for ', length(unlist(snp.blocks)), ' SNPs in ', nsnpblock, ' blocks...')

    beta <- foreach(k = 1:nsnpblock, .combine = rbind, .inorder = FALSE, .multicombine = TRUE) %do% {
        if(verbose) message('    Running block ', k, '...')
        # read genotype data for the block
        G <- .readGeno(gdsobj, sample.include, snp.index = snp.blocks[[k]])

        # calculate ISAF betas
        .calcISAFBeta(G = G, VVtVi = VVtVi)
    }
    ### rather than returning and rbinding here, we may want to be writing the output to something

    return(beta)
}


### exported function that does the pcrelate estimation for a (pair of) sample block(s)
pcrelateSampBlock <- function(gdsobj, betaobj, pcs, sample.include.block1, sample.include.block2,
                              scale = c('overall', 'variant', 'none'), ibd.probs = TRUE,
                              maf.thresh = 0.01, maf.bound.method = c('filter', 'truncate'),
                              verbose = TRUE){

    scale <- match.arg(scale)
    maf.bound.method <- match.arg(maf.bound.method)
    
    # create (joint) PC matrix and indices
    sample.include <- unique(c(sample.include.block1, sample.include.block2))
    V <- .createPCMatrix(pcs = pcs, sample.include = sample.include)
    idx <- which(rownames(V) %in% sample.include.block1)
    jdx <- which(rownames(V) %in% sample.include.block2)
    oneblock <- setequal(idx, jdx)
    ### slight inefficiency above because we create V for samples in block1 for each block2 when we don't have to if we are running serially; 
    ### but this seems more straightforward to parallelize

    snp.blocks <- .snpBlocks(gdsobj)
    nsnpblock <- length(snp.blocks)

    if(verbose) message('Running PC-Relate analysis using ', length(unlist(snp.blocks)), ' SNPs in ', nsnpblock, ' blocks...')
    # compute estimates for each variant block; sum along the way
    matList <- foreach(k = 1:nsnpblock, .combine = .matListCombine, .inorder = FALSE, .multicombine = FALSE) %do% {
        if(verbose) message('    Running block ', k, '...')
        # read genotype data for the block
        G <- .readGeno(gdsobj, sample.include, snp.index = snp.blocks[[k]])

        # load betas for the current block of variants
        beta.block <- betaobj[colnames(G), , drop = FALSE]
        ### this line of code will probably be different if we save the betas; need to load correct betas

        # calculate PC-Relate estimates
        .pcrelateVarBlock(	G = G, beta = beta.block, V = V, idx = idx, jdx = jdx, scale = scale, ibd.probs = ibd.probs, maf.thresh = maf.thresh, maf.bound.method = maf.bound.method)
    }

    # take ratios to compute final estimates
    estList <- .pcrelateCalcRatio(matList = matList, scale = scale, ibd.probs = ibd.probs)
    rm(matList)
    
    # cast to data.tables
    if(oneblock){
        # self table
        kinSelf <- data.table(ID = rownames(estList$kin), f = 2*diag(estList$kin) - 1, nsnp = diag(estList$nsnp))
        setkeyv(kinSelf, 'ID')
        # between samples table
        kinBtwn <- .estListToDT(estList, drop.lower = TRUE)
    }else{
        kinSelf <- NULL
        # between samples table
        kinBtwn <- .estListToDT(estList, drop.lower = FALSE)
    }
    rm(estList)
    setkeyv(kinBtwn, c('ID1', 'ID2'))
    
    return(list(kinSelf = kinSelf, kinBtwn = kinBtwn))
}
