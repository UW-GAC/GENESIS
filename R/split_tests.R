
#' Title Generate ID list for use in split tests
#'
#' @param data dataframe containing IDs matching genetic data sample IDs, and variable to divide groups by
#' @param id_var name of variable containing IDs
#' @param group_var name of variable to divide groups by
#'
#' @return a named list of individuals separated by group
#' @export
#'
#' @examples
gen_id_list <- function(data, id_var, group_var){
  if (!inherits(data, "data.frame")) stop("Input data is not a dataframe.") 
  groups <- unique(as.character(data[[group_var]]))
  id_list <- vector(mode="list", length=length(groups))
  names(id_list) <- groups
  for (g in groups){
    id_list[[g]] <- data[[id_var]][data[[group_var]]==g]
  }
  id_list
}



nullModelSplit <- function(nullmod, id_list, keep_all=TRUE){
  id_length <- ifelse(keep_all==TRUE, length(id_list)+1, length(id_list))
  split_nullmod <- vector(mode="list", length=id_length)
  list_names <- names(id_list)
  if (keep_all==TRUE) list_names <- c('all', names(id_list))
  names(split_nullmod) <- list_names
  keep_obj <- c('family', 'sample.id',  'outcome', 'workingY', 'resid.conditional', "resid.marginal", 'fitted.values') #"resid.marginal",  not needed?
  for (obj in names(nullmod)){
    if (!obj %in% keep_obj){
      nullmod[[obj]] <- NULL
    }
  }
  
  if (keep_all) split_nullmod[['all']] <- nullmod
  for (i in 1:length(id_list)){
    split_index <- ifelse(keep_all==TRUE, i+1, i)
    ids <- id_list[[i]]
    nullmod_subset <- nullmod
    idx <- nullmod$sample.id %in% ids
    for (v in c('sample.id',  'outcome', 'workingY', 'resid.conditional', 'fitted.values', "workingY")) {
      nullmod_subset[[v]] <- nullmod[[v]][idx]
    }
    split_nullmod[[split_index]] <- nullmod_subset
  }
  split_nullmod
}


setGeneric("assocTestSingle_split", function(gdsobj, ...) standardGeneric("assocTestSingle_split"))

setMethod("assocTestSingle_split",
          "SeqVarIterator",
  function(gdsobj, null.model, id_list, test=c("BinomiRare", "CMP"), 
          sparse=TRUE, imputed=FALSE, male.diploid=TRUE, genome.build=c("hg19", "hg38"),
          max.alt.freq=NULL, keep_all=TRUE, verbose=TRUE
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
  null.model_list <- nullModelSplit(null.model, id_list, keep_all=keep_all)
  
  
  # initialize results objects
  all_res <- rep( list(vector("list", length=n.iter)), length(null.model_list) ) #vector(mode = "list", length = length(id_list))
  names(all_res) <- names(null.model_list)
  
  
  #  geno_list <- vector(mode="list", length=length(id_list)) ###doing this is incorrect due to iterator behavior
  group_index_list <- vector(mode="list", length=length(null.model_list))
  
  
  for (i in 1:length(null.model_list)){
    current_ids <- null.model_list[[i]][["sample.id"]]
    cur_group_index <-match(current_ids, seqGetData(gdsobj, "sample.id"))
    group_index_list[[i]] <- cur_group_index
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
    
    for (grp.ind in 1:length(null.model_list)){
      cur_group_index <- group_index_list[[grp.ind]]
      current_nullmod <-null.model_list[[grp.ind]]
      message('length of current group ids: ', length(null.model_list[[grp.ind]][["sample.id"]]))
      current_geno <- geno[cur_group_index,,drop=FALSE]
      #     sample.index.grp <- which(is.element(rownames(current_geno), cur_group_ids)) ## is this a correct way to index these?? current_geo doesn't have rownames
      
      # allele frequency
      freq <- .alleleFreq(gdsobj, geno, sample.index=sample.index,
                          male.diploid=male.diploid, genome.build=genome.build)
      
      # take note of number of non-missing samples
      n.obs <- colSums(!is.na(current_geno))
      # filter monomorphic variants
      keep <- .filterMonomorphic(current_geno, count=n.obs, freq=freq$freq, imputed=imputed)
      if (!is.null(max.alt.freq)){
        keep <- keep & (freq$freq <= max.alt.freq)
      }
      
      if (!all(keep)) {
        current.var.info <- var.info[keep,,drop=FALSE] 
        current_geno <- current_geno[,keep,drop=FALSE]
        n.obs <- n.obs[keep]
        freq <- freq[keep,,drop=FALSE]
      } else {
        current.var.info <- var.info
      }
      
      # mean impute missing values
      if (any(n.obs < nrow(current_geno))) {
        current_geno <- .meanImpute(current_geno, freq$freq)
      }
      
      if (ncol(current_geno)==0){
        all_res[[grp.ind]][[i]]<- NULL
      } else{
      # do the test
      assoc <- testGenoSingleVar(current_nullmod, G=current_geno, test=test, calc_score=FALSE)
      # set monomorphs to NA - do we want to skip testing these to save time? ###not sure why this is needed when monomorphics were already filtered
      assoc[freq %in% c(0,1),] <- NA
      all_res[[grp.ind]][[i]]<- cbind(current.var.info, n.obs, freq, assoc)
      }
    }
    
    if (verbose & n.iter > 1 & i %% set.messages == 0) {
      message(paste("Iteration", i , "of", n.iter, "completed"))
    }
    i <- i + 1
    iterate <- iterateFilter(gdsobj, verbose=FALSE)
  }
  
  for (grp.ind in 1:length(all_res)){
    all_res[[grp.ind]] <- do.call(rbind, all_res[[grp.ind]])
  }
  all_res
}
)

save_split_results <- function(res_list, output_prefix=NULL){
    if (!is.null(output_prefix)){
      for (name in names(res_list)){
        output_path <- paste0(output_prefix, '_', name, '.RData')
        save(res_list[[name]], file=output_path)
      }
  } else {
      output_path <- paste0(name, '_results.RData')
      save(res_list[[name]], file=output_path)
  }
}


setGeneric("assocTestAggregate_split", function(gdsobj, ...) standardGeneric("assocTestAggregate_split"))

setMethod("assocTestAggregate_split",
          "SeqVarIterator",
          function(gdsobj, null.model, id_list, AF.max=1,
          #         weight.beta=c(1,1), ##all weights are set to 1
                   burden.test=c("BinomiRare", "CMP"), keep_all=TRUE,
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
            null.model_list <- nullModelSplit(null.model, id_list, keep_all=keep_all)
            n.iter <- length(variantFilter(gdsobj))
            # initialize results objects
            all_res <- rep( list(vector("list", length=n.iter)), length(null.model_list) ) #vector(mode = "list", length = length(id_list))
            names(all_res) <- names(null.model_list)
            
            all_res.var <- rep( list(vector("list", length=n.iter)), length(null.model_list) )
            names(all_res.var) <- names(null.model_list)
            
            #  geno_list <- vector(mode="list", length=length(id_list)) ###doing this is incorrect due to iterator behavior
            group_index_list <- vector(mode="list", length=length(null.model_list))
            
            
            for (i in 1:length(null.model_list)){
              current_ids <- null.model_list[[i]][["sample.id"]]
              cur_group_index <-match(current_ids, seqGetData(gdsobj, "sample.id"))
              group_index_list[[i]] <- cur_group_index
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
              for (grp.ind in 1:length(null.model_list)){
                cur_group_index <- group_index_list[[grp.ind]]
                current_nullmod <-null.model_list[[grp.ind]]
                message('length of current group ids: ', length(null.model_list[[grp.ind]][["sample.id"]]))
                current_geno <- geno[cur_group_index,,drop=FALSE]
                
                n.obs <- colSums(!is.na(current_geno))
                ###uncertain if this calculates correctly. initial code: freq <- .alleleFreq(gdsobj, geno, variant.index=index, sample.index=sample.index)
                freq <- .alleleFreq(gdsobj, current_geno, variant.index=index, sample.index=cur_group_index,
                  male.diploid=male.diploid, genome.build=genome.build)
                # number of non-missing samples
                
                # filter monomorphic variants
                keep <- .filterMonomorphic(current_geno, count=n.obs, freq=freq$freq, imputed=imputed)
                
                
                # exclude variants with freq > max
                keep <-  keep & freq$freq <= AF.max
                if (!all(keep)) {
                  current.var.info <- var.info[keep,,drop=FALSE]
                  current_geno <- current_geno[,keep,drop=FALSE]
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
                n.alt <- sum(current_geno, na.rm=TRUE)
                
                # number of samples with observed alternate alleles > 0
                n.sample.alt <- sum(rowSums(current_geno, na.rm=TRUE) > 0)
                
                all_res[[grp.ind]][[i]] <- data.frame(n.site, n.alt, n.sample.alt)
                all_res.var[[grp.ind]][[i]] <- cbind(current.var.info , n.obs, freq, weight)
                
                if (n.site > 0) {
                  # mean impute missing values
                  if (any(n.obs < nrow(current_geno))) {
                    current_geno <- .meanImpute(current_geno, freq$freq)
                  }
                  
                  # do the test
                  assoc <- testVariantSet(current_nullmod, G=current_geno, weights=weight, test="Burden", 
                                          burden.test=burden.test)
                  all_res[[grp.ind]][[i]] <- cbind(all_res[[grp.ind]][[i]], assoc)
                }
              }
              if (verbose & n.iter > 1 & i %% set.messages == 0) {
                message(paste("Iteration", i , "of", n.iter, "completed"))
              }
              i <- i + 1
              iterate <- iterateFilter(gdsobj, verbose=FALSE)
            }

            for (grp.ind in 1:length(all_res)){
              all_res[[grp.ind]] <- list(results=bind_rows(all_res[[grp.ind]]), variantInfo=all_res.var[[grp.ind]])
              all_res[[grp.ind]] <- .annotateAssoc(gdsobj, all_res[[grp.ind]])
            }
            all_res
          })


