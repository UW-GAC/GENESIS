pcair <- function(gdsobj,
                   kinobj,
                   divobj,
                   kin.thresh = 2^(-11/2),
                   div.thresh = -2^(-11/2),
                   unrel.set = NULL,
                   sample.include = NULL,
                   snp.include = NULL,
                   num.thread = 1L,
                   verbose = TRUE,
                   ...) {

    part <- pcairPartition(kinobj = kinobj, divobj = divobj,
                           kin.thresh = kin.thresh, div.thresh = div.thresh,
                           unrel.set = unrel.set, sample.include = sample.include,
                           verbose = verbose)

    if(verbose){
        message(paste("Unrelated Set:", length(part$unrels), "Samples",
                      "\nRelated Set:", length(part$rels), "Samples"))
    }

    if(verbose){  message("Performing PCA on the Unrelated Set...")  }
    pca.unrel <- snpgdsPCA(gdsobj, sample.id = part$unrels, snp.id = snp.include,
                           num.thread = num.thread, verbose = verbose, ...)

    if(length(part$rels > 0)){
        method <- "PC-AiR"
        if(verbose){  message("Predicting PC Values for the Related Set...") }
        snp.load <- snpgdsPCASNPLoading(pca.unrel, gdsobj = gdsobj,
                                        num.thread = num.thread, verbose = verbose)
        samp.load <- snpgdsPCASampLoading(snp.load, gdsobj = gdsobj,
                                          sample.id = part$rels,
                                          num.thread = num.thread, verbose = verbose)

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
    }
    eigenval <- pca.unrel$eigenval[1:ncol(pca.unrel$eigenvect)]
    
    # output object
    out <- list(vectors = eigenvect,
                values = eigenval,
                rels = part$rels,
                unrels = part$unrels,
                kin.thresh = kin.thresh,
                div.thresh = -abs(div.thresh),
                nsamp = length(unlist(part)),
                nsnps = length(pca.unrel$snp.id),
                varprop = pca.unrel$varprop,
                call = match.call(),
                method = method)
    class(out) <- "pcair"

    return(out)
}
