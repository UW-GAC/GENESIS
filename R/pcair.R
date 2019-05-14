setGeneric("pcair", function(gdsobj, ...) standardGeneric("pcair"))

setOldClass("gds.class")
setMethod("pcair",
          "gds.class",
          function(gdsobj, 
                   kinobj = NULL,
                   divobj = NULL,
                   kin.thresh = 2^(-11/2),
                   div.thresh = -2^(-11/2),
                   unrel.set = NULL,
                   sample.include = NULL,
                   snp.include = NULL,
                   num.cores = 1L,
                   verbose = TRUE,
                   ...) {
              .pcair(gdsobj = gdsobj, 
                     kinobj = kinobj,
                     divobj = divobj,
                     kin.thresh = kin.thresh,
                     div.thresh = div.thresh,
                     unrel.set = unrel.set,
                     sample.include = sample.include,
                     snp.include = snp.include,
                     num.cores = num.cores,
                     verbose = verbose,
                     ...)
          })

setOldClass("SNPGDSFileClass")
setMethod("pcair",
          "SNPGDSFileClass",
          function(gdsobj, ...) {
              .pcair(gdsobj, ...)
          })

setMethod("pcair",
          "GdsGenotypeReader",
          function(gdsobj, ...) {
              pcair(gdsobj@handler, ...)
          })

setMethod("pcair",
          "GenotypeData",
          function(gdsobj, ...) {
              pcair(gdsobj@data, ...)
          })

setMethod("pcair",
          "SeqVarGDSClass",
          function(gdsobj, ...) {
              filt <- seqGetFilter(gdsobj)
              out <- .pcair(gdsobj, ...)
              seqSetFilter(gdsobj,
                           sample.sel=filt$sample.sel,
                           variant.sel=filt$variant.sel,
                           verbose=FALSE)
              out
          })

.pcair <- function(gdsobj,
                   kinobj = NULL,
                   divobj = NULL,
                   kin.thresh = 2^(-11/2),
                   div.thresh = -2^(-11/2),
                   unrel.set = NULL,
                   sample.include = NULL,
                   snp.include = NULL,
                   num.cores = 1L,
                   verbose = TRUE,
                   ...) {

    if(!is.null(kinobj) & !is.null(divobj)){
      part <- pcairPartition(kinobj = kinobj, divobj = divobj,
                             kin.thresh = kin.thresh, div.thresh = div.thresh,
                             unrel.set = unrel.set, sample.include = sample.include,
                             verbose = verbose)

    }else if(is.null(kinobj) & is.null(divobj)){
      part <- .pcairPartitionUser(gdsobj = gdsobj, unrel.set = unrel.set, 
                                  sample.include = sample.include, verbose = verbose)

    }else{
      stop('kinobj and divobj must either both be specified, or both be NULL')
    }

    if(verbose){
        message(paste("Unrelated Set:", length(part$unrels), "Samples",
                      "\nRelated Set:", length(part$rels), "Samples"))
    }

    if(verbose){  message("Performing PCA on the Unrelated Set...")  }
    ## suppressMessages gets rid of the hint to use snpgdsOpen
    pca.unrel <- suppressMessages(
        snpgdsPCA(gdsobj, sample.id = part$unrels, snp.id = snp.include,
                  num.thread = num.cores, verbose = verbose, ...)
    )

    if(length(part$rels > 0)){
        method <- "PC-AiR"
        if(verbose){  message("Predicting PC Values for the Related Set...") }
        snp.load <- suppressMessages(
            snpgdsPCASNPLoading(pca.unrel, gdsobj = gdsobj,
                                num.thread = num.cores, verbose = verbose)
        )
        samp.load <- suppressMessages(
            snpgdsPCASampLoading(snp.load, gdsobj = gdsobj,
                                 sample.id = part$rels,
                                 num.thread = num.cores, verbose = verbose)
        )

        # combine unrelated and related PCs and order as in GDS file
        eigenvect <- rbind(pca.unrel$eigenvect, samp.load$eigenvect)
        rownames(eigenvect) <- c(as.character(pca.unrel$sample.id),
                                 as.character(samp.load$sample.id))
        sample.id <- as.character(read.gdsn(index.gdsn(gdsobj, "sample.id")))
        sample.id <- intersect(sample.id, rownames(eigenvect))
        samp.ord <- match(sample.id, rownames(eigenvect))
        eigenvect <- eigenvect[samp.ord,]

    } else {
        method <- "Standard PCA"
        eigenvect <- pca.unrel$eigenvect
        rownames(eigenvect) <- as.character(pca.unrel$sample.id)
    }
    eigenval <- pca.unrel$eigenval[1:ncol(pca.unrel$eigenvect)]
    
    # output object
    out <- list(vectors = eigenvect,
                values = eigenval,
                rels = part$rels,
                unrels = part$unrels,
                kin.thresh = kin.thresh,
                div.thresh = div.thresh,
                sample.id = rownames(eigenvect),
                nsamp = nrow(eigenvect),
                nsnps = length(pca.unrel$snp.id),
                varprop = pca.unrel$varprop,
                call = match.call(),
                method = method)
    class(out) <- "pcair"

    return(out)
}
