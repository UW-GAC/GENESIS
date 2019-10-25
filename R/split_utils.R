###current functions are only writen for assocTestSingleSplit results##
.detectVarNames <- function(dat, varname){
  if (class(dat)=="list"){
    any(unlist(lapply(dat, function(x){varname %in% names(x)}))==TRUE)
  } else if (inherits(dat, "data.frame")){
    any((varname %in% names(dat))==TRUE)
  }
}

matchSignifHits <- function(res.list, threshold, return.df=FALSE, variant.id.var=NULL){
  matched.results <- vector(mode="list", length=length(res.list))
  names(matched.results) <- names(res.list)
  subset.signif <- lapply(res.list, function(x) filter(x, pval < threshold))
  
  if (is.null(variant.id.var)){
    if (!.detectVarNames(res.list, c("variant.id", "variantID"))) stop("specify variant.id.var")
    if (.detectVarNames(res.list, "variant.id")) {
      variant.id.var <- "variant.id"
    } else if (.detectVarNames(res.list, "variantID")){
      variant.id.var <- "variantID"
    }
  }
  for (group in names(subset.signif)){
    if (nrow(subset.signif[[group]]) > 0){
      current.list <- res.list
      current.list[[group]] <- NULL
      subset.signif[[group]]$ref.group=group
      subset.signif[[group]]$signif.group=group
      matched.list <- lapply(current.list, function(x) filter(x, !!sym(variant.id.var) %in% subset.signif[[group]][[variant.id.var]]))
      for (g in names(matched.list)){
        if (nrow(matched.list[[g]]) > 0) {
          matched.list[[g]]$ref.group <- g
          matched.list[[g]]$signif.group <- group
        }
      } 
      matched.list <- bind_rows(matched.list)
      matched.results[[group]] <- bind_rows(subset.signif[[group]], matched.list) %>% arrange(!!sym(variant.id.var))
      
    } else {
      matched.results[[group]] <- NULL
    }
   
  }
  if (return.df) matched.results <- bind_rows(matched.results)
  matched.results
}


.filterIncomplete <- function(dat, total.groups, variant.id.var=NULL){
  if (is.null(variant.id.var)){
    if (!.detectVarNames(dat, c("variant.id", "variantID"))) stop("specify variant.id.var")
    if (.detectVarNames(dat, "variant.id")) {
      variant.id.var <- "variant.id"
    } else if (.detectVarNames(dat, "variantID")){
      variant.id.var <- "variantID"
    }
  }
  dat <- dat %>% filter(ref.group!="all") %>% group_by(!!sym(variant.id.var), signif.group) %>% 
    summarise(n=n()) %>% filter(n < total.groups) #, groups=glue::glue_collapse(unique(ref.group), sep=','), total_carriers=sum(n.carrier))
  dat[[variant.id.var]]
}


findIncompleteHits <- function(matched.results, n.groups=NULL, variant.id.var=NULL){
  if (class(matched.results)=="list"){
    if (is.null(n.groups)) n.groups <- sum(names(matched.results)!='all')
    incomplete.results <- lapply(matched.results, .filterIncomplete, total.groups=n.groups, variant.id.var=variant.id.var)
    incomplete.variants <- c()
    for (i in seq_along(incomplete.results)){
      incomplete.variants <- union(incomplete.variants, incomplete.results[[i]])
    }
  } else if (inherits(matched.results, "data.frame")){
    n.groups <- length(unique(matched.results$ref.group[matched.results$ref.group!="all"]))
    incomplete.results <- .filterIncomplete(matched.results, total.groups=n.groups, variant.id.var=variant.id.var)
    incomplete.variants <- unique(incomplete.results)
  }
  incomplete.variants
}



mergeNullModelBR <- function(nullmod.list, gdsfile){
  n.null <- length(nullmod.list)
  outcome.list <- fitted.values.list <- vector(mode = "list", length = n.null)
  
  for (i in 1:n.null){
    nullmod.i <- nullmod.list[[i]]
    
    if (nullmod.i$family$mixedmodel) { ## if this is a mixed model, used conditional probabilities
      phat <- expit(nullmod.i$workingY - nullmod.i$resid.conditional)    
    } else{ ## not a mixed model
      phat <- nullmod.i$fitted.values
    }
    names(phat) <- rownames(nullmod.i$model.matrix)
    fitted.values.list[[i]] <- phat
    
    outcome <- nullmod.i$outcome
    names(outcome) <- rownames(nullmod.i$model.matrix)
    outcome.list[[i]] <- outcome
    
  }
  
  fitted.values <- do.call(c, fitted.values.list)
  outcome <- do.call(c, outcome.list)
  
  ## re-order according to the order on the gds file: 
  if (is.character(gdsfile)){
    gds <- seqOpen(gdsfile)
  } else {
    gds <- gdsfile
  }
  sample.id <- seqGetData(gds, "sample.id")
  
  if (is.character(gdsfile)) seqClose(gds)

  ids.both <- intersect(sample.id, names(outcome))
  outcome <- outcome[match(ids.both, names(outcome))]
  fitted.values <- fitted.values[match(ids.both, names(fitted.values))]
  
  ##  set up the new (tricked) object. It needs to pass the checks for
  ## binomiRare: to have family = "binomial", not be a mixed model (for unified pull of)
  ## probability vector)
  ## fitted.values would be the probabilities; outcome the vectof of disease statuses. 
  new.nullmod <- list(family = list(family = "binomial", mixedmodel = FALSE), fitted.values = fitted.values, outcome = outcome, sample.id = ids.both)
  
  return(new.nullmod)
}




###need to subset and create iterator again

recreateIterator <- function(gds, annot, incomplete.variants, block.size=1024, verbose=TRUE){
  if (length(incomplete.variants)==0) stop("There are no variants to filter on. The gds object remains unchanged.")
  seqResetFilter(gds, verbose=verbose)
  seqSetFilter(gds, variant.id = incomplete.variants, verbose=verbose)
  seqData <- SeqVarData(gds, sampleData=annot)
  SeqVarBlockIterator(seqData, variantBlock=block.size, verbose=verbose)
}


runSplitSubset <- function(gdsobj, null.model, id.list, test=c("BinomiRare", "CMP"), remove.groups=NULL,
                   sparse=TRUE, imputed=FALSE, keep.all=TRUE, male.diploid=TRUE,  genome.build=c("hg19", "hg38"), verbose=TRUE) {
            message('running split version of assocTestSingle')
            test <- match.arg(test)
            if (!is.null(remove.groups)){
              for (g in remove.groups){
                  id.list[[g]] <- NULL
                }
            }
            test <- match.arg(test)
            # don't use sparse matrices for imputed dosages
            if (imputed) sparse <- FALSE
            # coerce null.model if necessary
            if (sparse) null.model <- .nullModelAsMatrix(null.model)
            
            # filter samples to match null model
            sample.index <- .setFilterNullModel(gdsobj, null.model, verbose=verbose) ##this subsets iterator sample.id to only those in null model
            
            # check ploidy
            if (SeqVarTools:::.ploidy(gdsobj) == 1) male.diploid <- FALSE
            
            n.iter <- length(variantFilter(gdsobj))
            set.messages <- ceiling(n.iter / 100) # max messages = 100
            
            
            ###split null model and genotypes by group###
            null.model.list <- nullModelSplit(null.model, id.list, keep.all=keep.all)
            
            
            # initialize results objects
            all.res <- rep( list(vector("list", length=n.iter)), length(null.model.list) ) 
            names(all.res) <- names(null.model.list)
            
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
#                message('length of current group ids: ', length(null.model.list[[grp.ind]][["sample.id"]]))
                current.geno <- geno[cur.group.index,,drop=FALSE]
                
                # allele frequency
                freq <- .alleleFreq(gdsobj, current.geno, sample.index=cur.group.index, male.diploid=male.diploid, genome.build=genome.build)
                # take note of number of non-missing samples
                n.obs <- colSums(!is.na(current.geno))

                if (any(n.obs < nrow(current.geno))) {
                  current.geno <- .meanImpute(current.geno, freq$freq)
                }

                  # do the test
                  assoc <- testGenoSingleVar(current.nullmod, G=current.geno, test=test)
                  assoc[freq$freq %in% c(0,1),] <- NA
                  all.res[[grp.ind]][[i]]<- cbind(var.info, n.obs, freq, assoc)
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


