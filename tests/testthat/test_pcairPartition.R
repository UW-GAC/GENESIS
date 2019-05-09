context("pcairPartition tests")

.testKinGds <- function(x, file) {
    gds <- gdsfmt::createfn.gds(file)
    gdsfmt::add.gdsn(gds, "sample.id", colnames(x))
    gdsfmt::add.gdsn(gds, "kinship", x)
    gdsfmt::closefn.gds(gds)
    gdsfmt::openfn.gds(file)
}

.cleanup <- function(x, file) {
    gdsfmt::closefn.gds(x)
    unlink(file)
}

test_that("name errors", {
    # load saved matrix of KING-robust estimates
    data("HapMap_ASW_MXL_KINGmat")

    newMat <- HapMap_ASW_MXL_KINGmat
    colnames(newMat) <- rownames(newMat) <- NULL
    expect_error(pcairPartition(kinobj = newMat, divobj = newMat, verbose=FALSE),
                 "colnames must be provided for kinobj and divobj")

    expect_warning(pcairPartition(kinobj = HapMap_ASW_MXL_KINGmat, divobj = HapMap_ASW_MXL_KINGmat, unrel.set = 1:100, verbose=FALSE),
                 "some samples in unrel.set are not in kinobj or divobj")
})

test_that("kinobj matrix, divobj gds", {
    data("HapMap_ASW_MXL_KINGmat")
    gdsfile <- tempfile()
    divobj <- .testKinGds(HapMap_ASW_MXL_KINGmat, gdsfile)

    mypart <- pcairPartition(kinobj = HapMap_ASW_MXL_KINGmat, divobj = divobj, verbose=FALSE)
    expect_equal(length(mypart$rels),76)
    expect_equal(length(mypart$unrels),97)

    .cleanup(divobj, gdsfile)
})

test_that("kinobj gds, divobj matrix", {
    data("HapMap_ASW_MXL_KINGmat")
    gdsfile <- tempfile()
    kinobj <- .testKinGds(HapMap_ASW_MXL_KINGmat, gdsfile)

    mypart <- pcairPartition(kinobj = kinobj, divobj = HapMap_ASW_MXL_KINGmat, verbose=FALSE)
    expect_equal(length(mypart$rels),76)
    expect_equal(length(mypart$unrels),97)

    .cleanup(kinobj, gdsfile)
})

test_that("matrix and gds give same results", {
    data("HapMap_ASW_MXL_KINGmat")
    gdsfile <- tempfile()
    kinobj <- .testKinGds(HapMap_ASW_MXL_KINGmat, gdsfile)
    
    mypart.mat <- pcairPartition(kinobj = HapMap_ASW_MXL_KINGmat, divobj = HapMap_ASW_MXL_KINGmat, verbose=FALSE)
    mypart.gds <- pcairPartition(kinobj = kinobj, divobj = kinobj, verbose=FALSE)
    expect_equal(mypart.mat, mypart.gds)

    .cleanup(kinobj, gdsfile)
})

test_that("kinobj and divobj both Matrix", {
    data("HapMap_ASW_MXL_KINGmat")
    Mat <- Matrix(HapMap_ASW_MXL_KINGmat)

    mypart <- pcairPartition(kinobj = Mat, divobj = Mat, verbose=FALSE)
    expect_equal(length(mypart$rels),76)
    expect_equal(length(mypart$unrels),97)
})

test_that("apply on Matrix", {
    x <- Matrix(matrix(1, nrow=10, ncol=20, dimnames=list(1:10,1:20)))
    MARGIN <- 1
    FUN <- sum
    selection <- list(1:5, 1:10)
    chk <- apply(x[selection[[1]], selection[[2]]], MARGIN = MARGIN, FUN = FUN)
    expect_equal(.apply(x, MARGIN, FUN, selection), chk)
    
    MARGIN <- 2
    chk <- apply(x[selection[[1]], selection[[2]]], MARGIN = MARGIN, FUN = FUN)
    expect_equal(.apply(x, MARGIN, FUN, selection), chk)
})

.apply.no.blocks <- function(x, MARGIN, FUN, selection) {
    x <- x[selection[[1]], selection[[2]]]
    ans <- list()
    if (MARGIN == 1) {
        for (i in 1:nrow(x)) {
            ans[[i]] <- FUN(x[i,])
        }
        names(ans) <- rownames(x)
    } else if (MARGIN == 2) {
        for (i in 1:ncol(x)) {
            ans[[i]] <- FUN(x[,i])
        }
        names(ans) <- colnames(x)
    } else {
        stop("MARGIN must be 1 or 2")
    }
    simplify2array(ans)
}

test_that("apply on big Matrix", {
    x <- Matrix(matrix(1, nrow=1000, ncol=2000, dimnames=list(1:1000,1:2000)))
    MARGIN <- 1
    FUN <- sum
    selection <- list(1:1000, 1:2000)
    chk <- apply(x[selection[[1]], selection[[2]]], MARGIN = MARGIN, FUN = FUN)
    expect_equal(.apply(x, MARGIN, FUN, selection, maxelem=5e5), chk)
    ## takes too long
    ## expect_equal(.apply(x, MARGIN, FUN, selection, maxelem=5e5),
    ##              .apply.no.blocks(x, MARGIN, FUN, selection))
    
    MARGIN <- 2
    chk <- apply(x[selection[[1]], selection[[2]]], MARGIN = MARGIN, FUN = FUN)
    expect_equal(.apply(x, MARGIN, FUN, selection, maxelem=5e5), chk)
})


test_that("no unrelated", {
    data("HapMap_ASW_MXL_KINGmat")
    Mat <- Matrix(HapMap_ASW_MXL_KINGmat)
    thresh <- min(Mat) - 0.1

    expect_error(pcairPartition(kinobj = Mat, divobj = Mat,
                                kin.thresh=thresh, div.thresh=thresh,
                                verbose=FALSE),
                 "All samples related")
})


test_that("no related", {
    Mat <- Diagonal(100)
    dimnames(Mat) <- list(1:100, 1:100)
    mypart <- pcairPartition(kinobj = Mat, divobj = Mat, verbose=FALSE)
    expect_true(is.null(mypart$rels))
    expect_equal(as.character(1:100), mypart$unrels)
})


test_that("all same number of relatives", {
    x <- matrix(c(0.5,0.25,0.25,0.5), nrow = 2)
    Mat <- bdiag(list(x,x,x,x,x))
    dimnames(Mat) <- list(1:10, 1:10)
    mypart <- pcairPartition(kinobj = Mat, divobj = Mat, verbose=FALSE)
    expect_equal(mypart$rels, as.character(seq(1,9,2)))
    expect_equal(mypart$unrels, as.character(seq(2,10,2)))
})
