context("test utils")

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
