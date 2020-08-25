setGeneric("effectAllele", function(gdsobj, ...) standardGeneric("effectAllele"))

setMethod("effectAllele",
          "SeqVarGDSClass",
          function(gdsobj, variant.id=NULL) {
              if (!is.null(variant.id)) {
                  var.info <- applyMethod(gdsobj, variantInfo, variant=variant.id,
                                          alleles=TRUE, expanded=TRUE)
              } else {
                  var.info <- variantInfo(gdsobj, alleles=TRUE, expanded=TRUE)
              }
              ## filt <- seqGetFilter(gdsobj)
              ## seqSetFilter(gdsobj, variant.id=variant.id)
              ## var.info <- variantInfo(gdsobj, alleles=TRUE, expanded=TRUE)
              ## seqSetFilter(gdsobj, variant.sel=filt$variant.sel)
              var.info <- var.info[,c("variant.id", "alt", "ref")]
              names(var.info)[2:3] <- c("effect.allele", "other.allele")
              var.info
          })

setMethod("effectAllele",
          "GenotypeData",
          function(gdsobj, variant.id=NULL) {
              A <- getAlleleA(gdsobj)
              B <- getAlleleB(gdsobj)
              if (is.null(A)) {
                  A <- getAlleleA(gdsobj@data)
                  B <- getAlleleB(gdsobj@data)
              }
              var.info <- data.frame(variant.id=getSnpID(gdsobj),
                                     effect.allele=A,
                                     other.allele=B,
                                     stringsAsFactors=FALSE)
              if (!is.null(variant.id)) {
                  var.info <- var.info[var.info$variant.id %in% variant.id,]
              }
              var.info
          })
