context("test utils")
library(SeqVarTools)

test_that("filterMonomorphic - discrete", {
    nsamp <- 100
    ref <- rep(0, nsamp)
    het <- rep(1, nsamp)
    alt <- rep(2, nsamp)
    set.seed(500); ref.miss <- ref; ref.miss[sample(nsamp, 5)] <- NA
    set.seed(501); het.miss <- het; het.miss[sample(nsamp, 5)] <- NA
    set.seed(502); alt.miss <- alt; alt.miss[sample(nsamp, 5)] <- NA
    set.seed(503); ok <- sample(rbinom(nsamp, 1, 0.2))
    set.seed(504); ok.miss <- ok; ok.miss[sample(nsamp, 5)] <- NA
    all.miss <- rep(NA, nsamp)
    geno <- cbind(ref,het,alt,ref.miss,het.miss,alt.miss,ok,ok.miss,all.miss)
    count <- colSums(!is.na(geno))
    freq <- 0.5*colMeans(geno, na.rm=TRUE)
    expect_equivalent(c(rep(FALSE, 6), rep(TRUE, 2), FALSE),
                      .filterMonomorphic(geno, count, freq))
})

test_that("filterMonomorphic - imputed", {
    nsamp <- 100
    set.seed(505); ref <- rep(0, nsamp) + runif(nsamp, 0, 1e-9)
    set.seed(506); het <- rep(1, nsamp) + runif(nsamp, -1e-9, 1e-9)
    set.seed(507); alt <- rep(2, nsamp) - runif(nsamp, 0, 1e-9)
    set.seed(508); ref.miss <- ref; ref.miss[sample(nsamp, 5)] <- NA
    set.seed(509); het.miss <- het; het.miss[sample(nsamp, 5)] <- NA
    set.seed(510); alt.miss <- alt; alt.miss[sample(nsamp, 5)] <- NA
    set.seed(511); ok <- runif(nsamp, 0, 2)
    set.seed(512); ok.miss <- ok; ok.miss[sample(nsamp, 5)] <- NA
    all.miss <- rep(NA, nsamp)
    geno <- cbind(ref,het,alt,ref.miss,het.miss,alt.miss,ok,ok.miss,all.miss)
    count <- colSums(!is.na(geno))
    freq <- 0.5*colMeans(geno, na.rm=TRUE)
    expect_equivalent(c(rep(FALSE, 6), rep(TRUE, 2), FALSE),
                      .filterMonomorphic(geno, count, freq, imputed=TRUE))
})

test_that("alleleFreq - autosomes", {
    svd <- .testData()
    freq <- alleleFrequency(svd)
    geno <- refDosage(svd)
    expect_equal(.alleleFreq(svd, geno)$freq, freq)
    seqClose(svd)
})

test_that("MAC - autosomes", {
    svd <- .testData()
    mac <- minorAlleleCount(svd)
    geno <- refDosage(svd)
    expect_equal(.alleleFreq(svd, geno)$MAC, mac)
    seqClose(svd)
})

test_that("alleleFreq - nosex", {
    svd <- .testData()
    sampleData(svd)$sex <- NULL
    freq <- alleleFrequency(svd)
    mac <- minorAlleleCount(svd)
    geno <- refDosage(svd)
    chk <- .alleleFreq(svd, geno)
    expect_equivalent(chk$freq, freq)
    expect_equivalent(chk$MAC, mac)
    seqClose(svd)
})

.testGdsXY <- function() {
    # make up file with sex chroms
    showfile.gds(closeall=TRUE, verbose=FALSE)
    gds.fn <- tempfile()
    invisible(file.copy(seqExampleFileName("gds"), gds.fn))
    gds <- openfn.gds(gds.fn, readonly=FALSE)
    node <- index.gdsn(gds, "chromosome")
    compression.gdsn(node, "")
    chr <- read.gdsn(node)
    chr[chr == 1] <- "X"
    chr[chr == 2] <- "Y"
    write.gdsn(node, chr)
    closefn.gds(gds)
    seqOptimize(gds.fn, target="chromosome", verbose=FALSE)
    gds <- seqOpen(gds.fn)
    sample.id <- seqGetData(gds, "sample.id")
    set.seed(56); sex <- sample(c("M","F"), replace=TRUE, length(sample.id))
    df <- data.frame(sample.id, sex, stringsAsFactors=FALSE)
    SeqVarData(gds, sampleData=Biobase::AnnotatedDataFrame(df))
}

.cleanupGds <- function(gds) {
    fn <- seqSummary(gds, check="none", verbose=FALSE)$filename
    seqClose(gds)
    unlink(fn)
}

## .test1KG_X <- function() {
##     gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
##     gdsfile <- system.file("extdata", "1KG_chrX.gds", package="SeqVarTools")
##     gds <- seqOpen(gdsfile)
##     data(sample_annotation_1KG)
##     SeqVarData(gds, sampleData=AnnotatedDataFrame(sample_annotation_1KG))
## }

test_that("alleleFreq - sex chrs", {
    svd <- .testGdsXY()
    freq <- alleleFrequency(svd)
    geno <- refDosage(svd)
    expect_equal(.alleleFreq(svd, geno)$freq, freq)
    .cleanupGds(svd)
})

test_that("MAC - sex chrs", {
    svd <- .testGdsXY()
    mac <- minorAlleleCount(svd)
    geno <- refDosage(svd)
    expect_equal(.alleleFreq(svd, geno)$MAC, round(mac))
    .cleanupGds(svd)
})

.test1KG_Y <- function() {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "1KG_chrY.gds", package="SeqVarTools")
    gds <- seqOpen(gdsfile)
    sample.id <- seqGetData(gds, "sample.id")
    df <- data.frame(sample.id, sex="M", stringsAsFactors=FALSE)
    svd <- SeqVarData(gds, sampleData=AnnotatedDataFrame(df))
    
    freq <- alleleFrequency(svd)
    geno <- refDosage(svd)
    chk <- .alleleFreq(svd, geno, male.diploid=FALSE)
    expect_equal(chk$freq, freq)
    
    mac <- minorAlleleCount(svd)
    expect_equal(chk$MAC, round(mac))
    seqClose(gds)
}

test_that("meanImpute", {
    n <- 1000
    #m <- 100000 takes too long
    m <- 1000
    geno <- matrix(rbinom(n*m, size = 2, prob = 0.1), nrow = n, ncol = m)
    
    miss <- sample(n*m, size = 0.1*n*m, replace = FALSE)
    geno[miss] <- NA
    
    freq <- 0.5*colMeans(geno, na.rm = TRUE)
    
    # original function
    x <- .meanImputeFn(geno, freq)
    
    # new function with matrix
    y <- .meanImpute(geno, freq)
    expect_equal(x, y)
    
    # new function with Matrix (one block)
    Geno <- Matrix(geno)
    y <- .meanImpute(Geno, freq)
    expect_equivalent(x, as.matrix(y))
    
    # new function with Matrix (multiple blocks)
    #n*m/2^25 # 3 blocks if m=100000
    y <- .meanImpute(Geno, freq, maxelem = 4e5)
    expect_equivalent(x, as.matrix(y))
})
