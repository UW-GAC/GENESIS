mat2gds <- function(mat, gdsfile) {
    ## if we have a Matrix object, coerce
    if (!is.matrix(mat)) mat <- as.matrix(mat)
    gds <- createfn.gds(gdsfile)
    add.gdsn(gds, "sample.id", colnames(mat))
    add.gdsn(gds, "kinship", mat)
    closefn.gds(gds)
}


kin2gds <- function(ibdobj, gdsfile) {
    gds <- createfn.gds(gdsfile)
    add.gdsn(gds, "sample.id", ibdobj$sample.id)
    add.gdsn(gds, "kinship", ibdobj$kinship)
    closefn.gds(gds)
}
