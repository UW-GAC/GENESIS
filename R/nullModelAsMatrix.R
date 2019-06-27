.checkMatrixType <- function(covMatList) {
    isMatrix <- sapply(covMatList, is, "Matrix")
    if (any(isMatrix) & !all(isMatrix)) {
        for (i in which(!isMatrix)) {
            covMatList[[i]] <- Matrix(covMatList[[i]], sparse=FALSE)
        }
    }
    covMatList
}


.nullModelAsMatrix <- function(null.model) {
    for (m in c("cholSigmaInv", "Ytilde", "resid", "CX", "CXCXI")) {
        if (is.matrix(null.model[[m]])) {
            null.model[[m]] <- Matrix(null.model[[m]], sparse=FALSE)
        }
    }
    null.model
}


.genoAsMatrix <- function(nullmod, G) {
    if (is(nullmod$cholSigmaInv, "Matrix") & is.matrix(G)) {
        G <- Matrix(G)
    }
    G
}
