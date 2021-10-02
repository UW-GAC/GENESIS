# assocMetaSingle() or metaTestSingle() or metaAnalysisSingle() or metaSingle()

# metaSingle <- function()

# metaAggregate <- function()

setGeneric("metaPrepScores", function(gdsobj, ...) standardGeneric("metaPrepScores"))

setMethod("metaPrepScores",
          "SeqVarIterator",
          function(gdsobj, null.model, score.cov = TRUE,
                   geno.coding=c("additive", "dominant", "recessive"),
                   sparse=TRUE, imputed=FALSE, male.diploid=TRUE, genome.build=c("hg19", "hg38"),
                   BPPARAM=bpparam(), verbose=TRUE) {

                     # better way to handle these checks?
                     if(class(gdsobj) == 'SeqVarBlockIterator' && score.cov){
                       stop('either set score.cov = FALSE when using a SeqVarBlockIterator object,
                       or use a SeqVarWindowIterator or SeqVarRangeIterator object')
                     }

                     if(class(gdsobj) == 'SeqVarWindowIterator' && !score.cov){
                       stop('either set score.cov = TRUE when using a SeqVarWindowIterator object,
                       or use a SeqVarBlockIterator or SeqVarRangeIterator object')
                     }

                     # genotype coding
                     geno.coding <- match.arg(geno.coding)

                     # don't use sparse matrices for imputed dosages
                     if (imputed) sparse <- FALSE

                     # Convert old null model format if necessary.
                     null.model <- .updateNullModelFormat(null.model)

                     # coerce null.model if necessary
                     if (sparse) null.model <- .nullModelAsMatrix(null.model)

                     # filter samples to match null model
                     sample.index <- .setFilterNullModel(gdsobj, null.model, verbose=verbose)

                     # get sex for calculating allele freq
                     sex <- validateSex(gdsobj)[sample.index]

                     # check ploidy
                     if (SeqVarTools:::.ploidy(gdsobj) == 1) male.diploid <- FALSE

                     if(verbose) message('Using ', bpnworkers(BPPARAM), ' CPU cores')

                     i <- 1
                     ITER <- function() {
                         iterate <- TRUE
                         if (i > 1) {
                             iterate <- iterateFilter(gdsobj, verbose=FALSE)
                         }
                         i <<- i + 1
                         if (!iterate) {
                             return(NULL)
                         }

                         # variant info
                         var.info <- variantInfo(gdsobj, alleles=TRUE, expanded=TRUE)
                         # chr:pos:ref:alt ID for easy merging across studies
                         VARID <- paste(paste0("chr", var.info$chr), var.info$pos, var.info$ref, var.info$alt, sep = ':')
                         var.info <- cbind(VARID, var.info[,1:5])
                         names(var.info)[5:6] <- c("other.allele", "effect.allele")

                         # read in genotype data
                         if (!imputed) {
                             geno <- expandedAltDosage(gdsobj, use.names=FALSE, sparse=sparse)[sample.index,,drop=FALSE]
                         } else {
                             geno <- imputedDosage(gdsobj, use.names=FALSE)[sample.index,,drop=FALSE]
                         }

                         chr <- chromWithPAR(gdsobj, genome.build=genome.build)

                         return(list(var.info=var.info, geno=geno, chr=chr))
                     }

                     # worker function
                     res <- bpiterate(ITER, .metaPrepScores, BPPARAM=BPPARAM,
                                      sex=sex, null.model=null.model,
                                      geno.coding=geno.coding,
                                      sparse=sparse, imputed=imputed,
                                      male.diploid=male.diploid,
                                      score.cov=score.cov)
                     .stopOnError(res)

                     # prepare output
                     if(score.cov){
                       res1 <- lapply(res, function(x) x[[1]])
                       res2 <- lapply(res, function(x) x[[2]])
                       rm(res)

                       # group names
                       grnames <- .metaGroupNames(gdsobj)
                       # variant info
                       names(res1) <- grnames
                       res1 <- as.data.frame(rbindlist(res1, use.names=TRUE, fill=TRUE, idcol = 'group.id')
                       # scores.cov
                       names(res2) <- grnames
                       list(variants = res1, scores.cov = res2)

                     }else{
                       out <- rbindlist(res)
                       out <- out[!duplicated(out)]
                       as.data.frame(out)
                     }
                 })


.metaPrepScores <- function(x, sex, null.model, geno.coding, sparse, imputed, male.diploid, score.cov, ...) {
  # prep the geno data
  x <- .prepGenoBlock(x, AF.max=1, geno.coding=geno.coding, imputed=imputed,
                      sex=sex, male.diploid=male.diploid)
  var.info <- x$var.info
  n.obs <- x$n.obs
  freq <- x$freq
  geno <- x$geno
  rm(x)

  # mean impute missing values
  if (any(n.obs < nrow(geno))) {
      geno <- .meanImpute(geno, freq$freq)
  }

  if (ncol(geno) == 0){
    res.i <- NULL
  } else {
    # calc score
    G <- .genoAsMatrix(null.model, geno)
    score <- as.vector(crossprod(G, null.model$fit$resid.PY)) # G^T P Y

    Gtilde <- calcGtilde(null.model, G)

    if(score.cov){
      # calc score covariance matrix
      V <- crossprod(Gtilde) # G^T P G
      score.SE <- sqrt(diag(V))
      colnames(V) <- rownames(V) <- var.info$VARID

      # # melt the covariance matrix to a data.table
      # V <- meltMatrix(GPG, drop.lower = TRUE, drop.diag = FALSE)
      # setnames(GPG, c('ID1', 'ID2'), c('VARID1', 'VARID2'))

      # output
      res.i <- list(cbind(var.info, n.obs, freq, Score = score, Score.SE = score.SE), V)

    }else{
      # calc score SE only
      GPG <- colSums(Gtilde^2) # vector of G^T P G (for each variant)
      score.SE <- sqrt(GPG)

      # output
      res.i <- cbind(var.info, n.obs, freq, Score = score, Score.SE = score.SE)
    }
  }

  return(res.i)
}


setGeneric(".metaGroupNames", function(gdsobj) standardGeneric(".metaGroupNames"))

setMethod(".metaGroupNames",
          "SeqVarWindowIterator",
          function(gdsobj) {
            gr <- variantRanges(gdsobj)
            grnames <- paste(paste0("chr", as.character(GenomicRanges::seqnames(gr))),
                            BiocGenerics::start(gr), BiocGenerics::end(gr), sep = "_")
            grnames
          })

setMethod(".metaGroupNames",
          "SeqVarRangeIterator",
          function(gdsobj) {
            # get variant group IDs
            gr <- variantRanges(gdsobj)
            grnames <- names(gr)
            if(is.null(grnames)){
              grnames <- paste(paste0("chr", as.character(GenomicRanges::seqnames(gr))),
                               BiocGenerics::start(gr), BiocGenerics::end(gr), sep = "_")
            }
            grnames
          })

setMethod(".metaGroupNames",
          "SeqVarListIterator",
          function(gdsobj, res1, res2) {
              stop('not implemented')
          })



#
# setGeneric(".metaPrepScoresOutput", function(gdsobj, res1, res2) standardGeneric(".metaPrepScoresOutput"))
#
# setMethod(".metaPrepScoresOutput",
#           "SeqVarWindowIterator",
#           function(gdsobj, res1, res2) {
#             gr <- variantRanges(gdsobj)
#             grnames <- paste(paste0("chr", as.character(GenomicRanges::seqnames(gr))),
#                             BiocGenerics::start(gr), BiocGenerics::end(gr), sep = "_")
#
#             # variant info
#             # for(i in 1:length(res1)){
#             #     res1[[i]] <- cbind(group.id = grnames[i], res1[[i]])
#             # }
#             # res1 <- as.data.frame(rbindlist(res1, use.names=TRUE, fill=TRUE))
#
#             # variant info
#             names(res1) <- grnames
#             res1 <- as.data.frame(rbindlist(res1, use.names=TRUE, fill=TRUE, idcol = 'group.id')
#
#             # scores.cov
#             names(res2) <- grnames
#
#             list(variants = res1, scores.cov = res2)
#           })
#
# setMethod(".metaPrepScoresOutput",
#           "SeqVarRangeIterator",
#           function(gdsobj, res1, res2) {
#             # get variant group IDs
#             gr <- variantRanges(gdsobj)
#             grnames <- names(gr)
#             if(is.null(grnames)){
#               grnames <- paste(paste0("chr", as.character(GenomicRanges::seqnames(gr))),
#                                BiocGenerics::start(gr), BiocGenerics::end(gr), sep = "_")
#             }
#
#             # variant info
#             # for(i in 1:length(res1)){
#             #     res1[[i]] <- cbind(group.id = grnames[i], res1[[i]])
#             # }
#             # res1 <- as.data.frame(rbindlist(res1, use.names=TRUE, fill=TRUE))
#
#             # variant info
#             names(res1) <- grnames
#             res1 <- as.data.frame(rbindlist(res1, use.names=TRUE, fill=TRUE, idcol = 'group.id')
#
#             # scores.cov
#             names(res2) <- grnames
#
#             list(variants = res1, scores.cov = res2)
#           })
#
