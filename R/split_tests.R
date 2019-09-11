
genIDList <- function(data, id.var, group.var){
  if (!inherits(data, "data.frame")) stop("Input data is not a dataframe.")
  if (!(id.var %in% names(data))) stop("The id.var provided is not a variable in data")
  if (!(group.var %in% names(data))) stop("The group.var provided is not a variable in data")
  data <- na.exclude(data[ , c(id.var, group.var)])
  groups <- unique(as.character(data[[group.var]]))
  id.list <- vector(mode="list", length=length(groups))
  names(id.list) <- groups
  for (g in groups){
    id.list[[g]] <- data[[id.var]][data[[group.var]]==g]
  }
  id.list
}



nullModelSplit <- function(nullmod, id.list, keep.all=TRUE){
  id.length <- ifelse(keep.all==TRUE, length(id.list)+1, length(id.list))
  split.nullmod <- vector(mode="list", length=id.length)
  list.names <- names(id.list)
  if (keep.all==TRUE) list.names <- c('all', names(id.list))
  names(split.nullmod) <- list.names
  keep.obj <- c('family', 'sample.id',  'outcome', 'workingY', 'resid.conditional', "resid.marginal", 'fitted.values') #"resid.marginal",  not needed?
  for (obj in names(nullmod)){
    if (!obj %in% keep.obj){
      nullmod[[obj]] <- NULL
    }
  }
  
  if (keep.all) split.nullmod[['all']] <- nullmod
  for (i in 1:length(id.list)){
    split.index <- ifelse(keep.all==TRUE, i+1, i)
    ids <- id.list[[i]]
    nullmod.subset <- nullmod
    idx <- nullmod$sample.id %in% ids
    for (v in c('sample.id',  'outcome', 'workingY', 'resid.conditional', 'fitted.values', "workingY")) {
      nullmod.subset[[v]] <- nullmod[[v]][idx]
    }
    split.nullmod[[split.index]] <- nullmod.subset
  }
  split.nullmod
}


setGeneric("assocTestSingleSplit", function(gdsobj, ...) standardGeneric("assocTestSingleSplit"))

setMethod("assocTestSingleSplit",
          "SeqVarIterator",
  function(gdsobj, null.model, id.list, test=c("BinomiRare", "CMP"), 
          sparse=TRUE, imputed=FALSE, male.diploid=TRUE, genome.build=c("hg19", "hg38"),
          AF.max=NULL, keep.all=TRUE, verbose=TRUE
                                  ) {
  message('running split version of assocTestSingle')
  test <- match.arg(test)
  # don't use sparse matrices for imputed dosages
  if (imputed) sparse <- FALSE
  # coerce null.model if necessary
  ### need to add options here to check null model  if (inherits(data, 'list')){ ??
  if (sparse) null.model <- .nullModelAsMatrix(null.model)
  
  
  # filter samples to match null model
  sample.index <- .setFilterNullModel(gdsobj, null.model, verbose=verbose) ##this subsets iterator sample.id to only those in null model
  
  
  if (SeqVarTools:::.ploidy(gdsobj) == 1) male.diploid <- FALSE
  
  n.iter <- length(variantFilter(gdsobj))
  set.messages <- ceiling(n.iter / 100) # max messages = 100
  
  
  ###split null model and genotypes by group###
  null.model.list <- nullModelSplit(null.model, id.list, keep.all=keep.all)
  
  
  # initialize results objects
  all.res <- rep( list(vector("list", length=n.iter)), length(null.model.list) ) #vector(mode = "list", length = length(id.list))
  names(all.res) <- names(null.model.list)
  
  
  #  geno.list <- vector(mode="list", length=length(id.list)) ###doing this is incorrect due to iterator behavior
  group.index.list <- vector(mode="list", length=length(null.model.list))
  
  
  for (i in 1:length(null.model.list)){
    current.ids <- null.model.list[[i]][["sample.id"]]
    cur.group.index <-match(current.ids, seqGetData(gdsobj, "sample.id"))
    group.index.list[[i]] <- cur.group.index
  }
  
  i <- 1
  iterate <- TRUE
  while (iterate) {
    if (!imputed) {
      geno <- expandedAltDosage(gdsobj, use.names=FALSE, sparse=sparse)[sample.index,,drop=FALSE]
    } else {
      geno <- imputedDosage(gdsobj, use.names=FALSE)[sample.index,,drop=FALSE]
    }
    var.info <- variantInfo(gdsobj, alleles=FALSE, expanded=TRUE)
    
    for (grp.ind in 1:length(null.model.list)){
      cur.group.index <- group.index.list[[grp.ind]]
      current.nullmod <-null.model.list[[grp.ind]]
      message('length of current group ids: ', length(null.model.list[[grp.ind]][["sample.id"]]))
      current.geno <- geno[cur.group.index,,drop=FALSE]
      #     sample.index.grp <- which(is.element(rownames(current.geno), cur.group.ids)) ## is this a correct way to index these?? current.geo doesn't have rownames
      
      # allele frequency
      freq <- .alleleFreq(gdsobj, current.geno, sample.index=sample.index,
                          male.diploid=male.diploid, genome.build=genome.build)
      
      # take note of number of non-missing samples
      n.obs <- colSums(!is.na(current.geno))
      # filter monomorphic variants
      keep <- .filterMonomorphic(current.geno, count=n.obs, freq=freq$freq, imputed=imputed)
      if (!is.null(AF.max)){
        keep <- keep & (freq$freq <= AF.max)
      }
      
      if (!all(keep)) {
        current.var.info <- var.info[keep,,drop=FALSE] 
        current.geno <- current.geno[,keep,drop=FALSE]
        n.obs <- n.obs[keep]
        freq <- freq[keep,,drop=FALSE]
      } else {
        current.var.info <- var.info
      }
      
      # mean impute missing values
      if (any(n.obs < nrow(current.geno))) {
        current.geno <- .meanImpute(current.geno, freq$freq)
      }
      
      if (ncol(current.geno)==0){
        all.res[[grp.ind]][[i]]<- NULL
      } else{
      # do the test
      assoc <- testGenoSingleVar(current.nullmod, G=current.geno, test=test, calc.score=FALSE)
      # set monomorphs to NA - do we want to skip testing these to save time? ###not sure why this is needed when monomorphics were already filtered
      assoc[freq %in% c(0,1),] <- NA
      all.res[[grp.ind]][[i]]<- cbind(current.var.info, n.obs, freq, assoc)
      }
    }
    
    if (verbose & n.iter > 1 & i %% set.messages == 0) {
      message(paste("Iteration", i , "of", n.iter, "completed"))
    }
    i <- i + 1
    iterate <- iterateFilter(gdsobj, verbose=FALSE)
  }
  
  for (grp.ind in 1:length(all.res)){
    all.res[[grp.ind]] <- do.call(rbind, all.res[[grp.ind]])
  }
  all.res
}
)

saveSplitResults <- function(res.list, output.prefix=NULL){
  for (name in names(res.list)){
    current.dat <- res.list[[name]]
    if (!is.null(output.prefix)){
        output.path <- paste0(output.prefix, '.', name, '.RData')
      } else {
      output.path <- paste0(name, '.results.RData')
    }
  save(current.dat, file=output.path)
  }
}

setGeneric("assocTestAggregateSplit", function(gdsobj, ...) standardGeneric("assocTestAggregateSplit"))

setMethod("assocTestAggregateSplit",
          "SeqVarIterator",
          function(gdsobj, null.model, id.list, AF.max=1,
          #         weight.beta=c(1,1), ##all weights are set to 1
                   burden.test=c("BinomiRare", "CMP"), keep.all=TRUE,
                   sparse=TRUE, imputed=FALSE, 
                  male.diploid=TRUE, genome.build=c("hg19", "hg38"),
                  verbose=TRUE) {
            burden.test <- match.arg(burden.test)
            # don't use sparse matrices for imputed dosages
            if (imputed) sparse <- FALSE
            
            # coerce null.model if necessary
            if (sparse) null.model <-  .nullModelAsMatrix(null.model)
            
            # filter samples to match null model
            sample.index <-  .setFilterNullModel(gdsobj, null.model, verbose=verbose)
            
            # do we need to match on alleles?
            match.alleles <- any(c("ref", "alt") %in% names(S4Vectors:::mcols(currentRanges(gdsobj)))) ##S4Vectors:::mcols
            
            # check ploidy
            if (SeqVarTools:::.ploidy(gdsobj) == 1) male.diploid <- FALSE
            
            ###split null model by group###
            null.model.list <- nullModelSplit(null.model, id.list, keep.all=keep.all)
            n.iter <- length(variantFilter(gdsobj))
            # initialize results objects
            all.res <- rep( list(vector("list", length=n.iter)), length(null.model.list) ) #vector(mode = "list", length = length(id.list))
            names(all.res) <- names(null.model.list)
            
            all.res.var <- rep( list(vector("list", length=n.iter)), length(null.model.list) )
            names(all.res.var) <- names(null.model.list)
            
            #  geno.list <- vector(mode="list", length=length(id.list)) ###doing this is incorrect due to iterator behavior
            group.index.list <- vector(mode="list", length=length(null.model.list))
            
            
            for (i in 1:length(null.model.list)){
              current.ids <- null.model.list[[i]][["sample.id"]]
              cur.group.index <-match(current.ids, seqGetData(gdsobj, "sample.id"))
              group.index.list[[i]] <- cur.group.index
            }
          
            i <- 1
            
            set.messages <- ceiling(n.iter / 100) # max messages = 100
            iterate <- TRUE
            while (iterate) {
              var.info <- variantInfo(gdsobj, alleles=match.alleles, expanded=TRUE)
              
              if (!imputed) {
                geno <- expandedAltDosage(gdsobj, use.names=FALSE, sparse=sparse)[sample.index,,drop=FALSE]
              } else {
                geno <- imputedDosage(gdsobj, use.names=FALSE)[sample.index,,drop=FALSE]
              }
              
              if (match.alleles) {
                index <-  .matchAlleles(gdsobj, var.info)
                var.info <- var.info[index,,drop=FALSE]
                geno <- geno[,index,drop=FALSE]
              } else {
                index <- NULL
              }
              
              ####start null model loop from here####
              for (grp.ind in 1:length(null.model.list)){
                cur.group.index <- group.index.list[[grp.ind]]
                current.nullmod <-null.model.list[[grp.ind]]
                message('length of current group ids: ', length(null.model.list[[grp.ind]][["sample.id"]]))
                current.geno <- geno[cur.group.index,,drop=FALSE]
                
                n.obs <- colSums(!is.na(current.geno))
                ###uncertain if this calculates correctly. initial code: freq <- .alleleFreq(gdsobj, geno, variant.index=index, sample.index=sample.index)
                freq <- .alleleFreq(gdsobj, current.geno, variant.index=index, sample.index=cur.group.index,
                  male.diploid=male.diploid, genome.build=genome.build)
                # number of non-missing samples
                
                # filter monomorphic variants
                keep <- .filterMonomorphic(current.geno, count=n.obs, freq=freq$freq, imputed=imputed)
                
                
                # exclude variants with freq > max
                keep <-  keep & freq$freq <= AF.max
                if (!all(keep)) {
                  current.var.info <- var.info[keep,,drop=FALSE]
                  current.geno <- current.geno[,keep,drop=FALSE]
                  n.obs <- n.obs[keep]
                  freq <- freq[keep,,drop=FALSE]
                } else {
                  current.var.info <- var.info
                }
                
                # weight
                weight <- seq(from=1, to=1, length.out=length(freq$freq))
                
                # number of variant sites
                n.site <- length(unique(current.var.info$variant.id))
                
                # number of alternate alleles
                n.alt <- sum(current.geno, na.rm=TRUE)
                
                # number of samples with observed alternate alleles > 0
                n.sample.alt <- sum(rowSums(current.geno, na.rm=TRUE) > 0)
                
                all.res[[grp.ind]][[i]] <- data.frame(n.site, n.alt, n.sample.alt)
                all.res.var[[grp.ind]][[i]] <- cbind(current.var.info , n.obs, freq, weight)
                
                if (n.site > 0) {
                  # mean impute missing values
                  if (any(n.obs < nrow(current.geno))) {
                    current.geno <- .meanImpute(current.geno, freq$freq)
                  }
                  
                  # do the test
                  assoc <- testVariantSet(current.nullmod, G=current.geno, weights=weight, test="Burden", 
                                          burden.test=burden.test)
                  all.res[[grp.ind]][[i]] <- cbind(all.res[[grp.ind]][[i]], assoc)
                }
              }
              if (verbose & n.iter > 1 & i %% set.messages == 0) {
                message(paste("Iteration", i , "of", n.iter, "completed"))
              }
              i <- i + 1
              iterate <- iterateFilter(gdsobj, verbose=FALSE)
            }

            for (grp.ind in 1:length(all.res)){
              all.res[[grp.ind]] <- list(results=bind_rows(all.res[[grp.ind]]), variantInfo=all.res.var[[grp.ind]])
              all.res[[grp.ind]] <- .annotateAssoc(gdsobj, all.res[[grp.ind]])
            }
            all.res
          })


