###current functions are only writen for assocTestSingle_split results##

match_signif_hits <- function(res_list, threshold, return_df=FALSE){
  matched_results <- vector(mode="list", length=length(res_list))
  names(matched_results) <- names(res_list)
  subset_signif <- lapply(res_list, function(x) filter(x, pval < threshold))
  for (group in names(subset_signif)){
    if (nrow(subset_signif[[group]]) > 0){
      current_list <- res_list
      current_list[[group]] <- NULL
      subset_signif[[group]]$ref_group=group
      subset_signif[[group]]$signif_group=group
      matched_list <- lapply(current_list, function(x) filter(x, variantID %in% subset_signif[[group]]$variantID))
      for (g in names(matched_list)){
        if (nrow(matched_list[[g]]) > 0) {
          matched_list[[g]]$ref_group <- g
          matched_list[[g]]$signif_group <- group
        }
      } 
      matched_list <- bind_rows(matched_list)
      matched_results[[group]] <- bind_rows(subset_signif[[group]], matched_list) %>% arrange(variantID)
      
    } else {
      matched_results[[group]] <- NULL
    }
   
  }
  if (return_df) matched_results <- bind_rows(matched_results)
  matched_results
}


filter_incomplete <- function(dat, total_groups){
  dat <- dat %>% filter(ref_group!="all") %>% group_by(variantID) %>% 
    summarise(n=n()) %>% filter(n < total_groups) #, groups=glue::glue_collapse(unique(ref_group), sep=','), total_carriers=sum(n.carrier))
  dat$variantID 
}


find_incomplete_hits <- function(matched_results, n_groups=NULL){
  if (class(matched_results)=="list"){
    if (is.null(n_groups)) n_groups <- length(matched_results)
    incomplete_results <- lapply(matched_results, filter_incomplete, total_groups=n_groups)
    incomplete_variants <- c()
    for (i in seq_along(incomplete_results)){
      incomplete_variants <- union(incomplete_variants, incomplete_results[[i]])
    }
  } else if (inherits(matched_results, "data.frame")){
    n_groups <- length(unique(matched_results$ref_group[matched_results$ref_group!="all"]))
    incomplete_results <- filter_incomplete(matched_results, total_groups=n_groups)
    incomplete_variants <- unique(incomplete_results)
  }
  incomplete_variants
}



###need to subset and create iterator again

recreate_iterator <- function(gds, annot, incomplete_variants, block.size=1024){
  seqResetFilter(gds)
  seqSetFilter(gds, variant.id = incomplete_variants)
  seqData <- SeqVarData(gds, sampleData=annot)
  SeqVarBlockIterator(seqData, variantBlock=block.size)
}


run_split_subset <- function(gdsobj, null.model, id_list, test=c("BinomiRare", "CMP"), remove_groups=NULL,
                   sparse=TRUE, imputed=FALSE, keep_all=TRUE, male.diploid=TRUE,  genome.build=c("hg19", "hg38"), verbose=TRUE) {
            message('running split version of assocTestSingle')
            test <- match.arg(test)
            if (!is.null(remove_groups)){
              for (g in remove_groups){
                  id_list[[g]] <- NULL
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
            null.model_list <- nullModelSplit(null.model, id_list, keep_all=keep_all)
            
            
            # initialize results objects
            all_res <- rep( list(vector("list", length=n.iter)), length(null.model_list) ) 
            names(all_res) <- names(null.model_list)
            
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
#                message('length of current group ids: ', length(null.model_list[[grp.ind]][["sample.id"]]))
                current_geno <- geno[cur_group_index,,drop=FALSE]
                
                # allele frequency
                freq <- .alleleFreq(gdsobj, current_geno, sample.index=cur_group_index, male.diploid=male.diploid, genome.build=genome.build)
                # take note of number of non-missing samples
                n.obs <- colSums(!is.na(current_geno))

                if (any(n.obs < nrow(current_geno))) {
                  current_geno <- .meanImpute(current_geno, freq$freq)
                }

                  # do the test
                  assoc <- testGenoSingleVar(current_nullmod, G=current_geno, test=test, calc_score=FALSE)
                  assoc[freq$freq %in% c(0,1),] <- NA
                  all_res[[grp.ind]][[i]]<- cbind(var.info, n.obs, freq, assoc)
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


