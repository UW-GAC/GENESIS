setGeneric(".apply", function(x, MARGIN, FUN, selection, ...) standardGeneric(".apply"))
setMethod(".apply",
          "matrix",
          function(x, MARGIN, FUN, selection) {
              apply(x[selection[[1]], selection[[2]]], MARGIN = MARGIN, FUN = FUN)
          })


setMethod(".apply",
          "gds.class",
          function(x, MARGIN, FUN, selection) {
              apply.gdsn(index.gdsn(x, 'kinship'), margin = MARGIN, FUN = FUN,
                         selection = selection)
          })


## apply coerces its argument with as.matrix, and this fails for a matrix
## with > 2^31 - 1 elements
setMethod(".apply",
          "Matrix",
          function (x, MARGIN, FUN, selection, maxelem = 2^30){
              
              # subset to selection
              x <- x[selection[[1]], selection[[2]]]
              
              # determine number of blocks needed
              nr <- as.numeric(nrow(x))
              nc <- as.numeric(ncol(x))
              nblock <- ceiling(nr*nc/maxelem)

              if(nblock > 1){
                  
                  if(MARGIN  == 1){
                      blocks <- unname(split(1:nr, cut(1:nr, nblock)))
                      ans <- lapply(blocks, function(b) {
                          # need to coerce output of apply to a list
                          as.list(apply(x[b,], 1, FUN))
                      })
                      
                  }else if(MARGIN == 2){
                      blocks <- unname(split(1:nc, cut(1:nc, nblock)))
                      ans <- lapply(blocks, function(b) {
                          # need to coerce output of apply to a list
                          as.list(apply(x[,b], 2, FUN))
                      })
                      
                  }else {
                      stop("MARGIN must be 1 or 2")
                  }
                  
                  # unlist the top level
                  ans <- unlist(ans, recursive = FALSE)
                  if (length(ans) == 0){
                      return(vector(mode(x[1,1]), length=0))
                  }
                  
                  # simplify further if possible
                  return(simplify2array(ans))
                  
              }else{
                  return(apply(x, MARGIN, FUN))
              }
          })


setGeneric(".countNonMissing", function(x, MARGIN) standardGeneric(".countNonMissing"))

setMethod(".countNonMissing",
          "matrix",
          function(x, MARGIN){
            if(MARGIN == 1){
              rowSums(!is.na(x))
            }else if(MARGIN == 2){
              colSums(!is.na(x))
            }
          })

setMethod(".countNonMissing",
          "Matrix",
          function(x, MARGIN){
            nr <- as.numeric(nrow(x))
            nc <- as.numeric(ncol(x))

            if(nr*nc < 2^31){
              if(MARGIN == 1){
                rowSums(!is.na(x))
              }else if(MARGIN == 2){
                colSums(!is.na(x))
              }
            }else{
              .apply(x, 
                    MARGIN = MARGIN, 
                    FUN = function(x){ sum(!is.na(x)) },
                    selection = list(1:nr, 1:nc))
            }
            })


setGeneric(".readSampleId", function(x) standardGeneric(".readSampleId"))
setMethod(".readSampleId",
          "matrix",
          function(x) {
              colnames(x)
          })

setMethod(".readSampleId",
          "Matrix",
          function(x) {
              colnames(x)
          })

setMethod(".readSampleId",
          "gds.class",
          function(x) {
              read.gdsn(index.gdsn(x, 'sample.id'))
          })

setMethod(".readSampleId",
          "GdsGenotypeReader",
          function(x) {
              getScanID(x)
          })

setMethod(".readSampleId",
          "GenotypeData",
          function(x) {
              getScanID(x)
          })

setMethod(".readSampleId",
          "SeqVarGDSClass",
          function(x) {
              seqGetData(x, "sample.id")
          })



setGeneric(".readGeno", function(gdsobj, ...) standardGeneric(".readGeno"))
setMethod(".readGeno",
          "SeqVarGDSClass",
          function(gdsobj, sample.include=NULL, snp.index=NULL, allele=c("alt", "ref")){
              seqSetFilter(gdsobj, sample.id=sample.include, variant.sel=snp.index,
                           verbose=FALSE)
              allele <- match.arg(allele)
              if (allele == "alt") {
                  return(altDosage(gdsobj))
              } else  {
                  return(refDosage(gdsobj))
              }
          })

setMethod(".readGeno",
          "GdsGenotypeReader",
          function(gdsobj, sample.include=NULL, snp.index=NULL){
              getGenotypeSelection(gdsobj, scanID=sample.include, snp=snp.index,
                                   transpose=TRUE, drop=FALSE)
          })

setMethod(".readGeno",
          "MatrixGenotypeReader",
          function(gdsobj, sample.include=NULL, snp.index=NULL){
              getGenotypeSelection(gdsobj, scanID=sample.include, snp=snp.index,
                                   transpose=TRUE, drop=FALSE)
          })

setMethod(".readGeno",
          "GenotypeData",
          function(gdsobj, ...){
              .readGeno(gdsobj@data, ...)
          })



setGeneric(".snpBlocks", function(gdsobj, ...) standardGeneric(".snpBlocks"))
setMethod(".snpBlocks",
          "SeqVarIterator",
          function(gdsobj) {
              variantFilter(gdsobj)
          })

setMethod(".snpBlocks",
          "GenotypeIterator",
          function(gdsobj) {
              snpFilter(gdsobj)
          })

