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
    expect_equal(.alleleFreq(svd, geno), freq)
    seqClose(svd)
})

test_that("MAC - autosomes", {
    svd <- .testData()
    mac <- minorAlleleCount(svd)
    geno <- refDosage(svd)
    expect_equal(.minorAlleleCount(svd, geno), mac)
    seqClose(svd)
})

.testGdsXY <- function() {
    # make up file with sex chroms
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
    set.seed(55); sex <- sample(c("M","F"), replace=TRUE, length(sample.id))
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
    expect_equal(.alleleFreq(svd, geno), freq)
    .cleanupGds(svd)
})

test_that("MAC - sex chrs", {
    svd <- .testGdsXY()
    mac <- minorAlleleCount(svd)
    geno <- refDosage(svd)
    expect_equal(.minorAlleleCount(svd, geno), round(mac))
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
    expect_equal(.alleleFreq(svd, geno, male.diploid=FALSE), freq)
    
    mac <- minorAlleleCount(svd)
    expect_equal(.minorAlleleCount(svd, geno, male.diploid=FALSE), round(mac))
    seqClose(gds)
}


