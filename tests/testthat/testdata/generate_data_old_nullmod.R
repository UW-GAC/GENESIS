# Run using an older version of GENESIS.

# Generate data.
set.seed(4)
n <- 20
dat <- data.frame(sample.id=letters[1:n],
                  outcome_linear=rnorm(n),
                  outcome_binary=sample(c(0, 1), n, replace = T),
                  b=c(rep("a", n/2), rep("b", n/2)),
                  group = rep(c("c", "d"), len.out = n),
                  stringsAsFactors=FALSE)
cov_mat <- crossprod(matrix(rnorm(n*n, sd=0.05), n, n))
dimnames(cov_mat) <- list(dat$sample.id, dat$sample.id)
dat <- AnnotatedDataFrame(dat)


# Linear
nm <- fitNullModel(dat, outcome="outcome_linear", covars="b", verbose=FALSE)
saveRDS(nm, file = "old_nullmod_lm.rds", version = 2)

# Linear mixed
nm <- fitNullModel(dat, outcome="outcome_linear", covars="b", verbose=FALSE, cov.mat = cov_mat)
saveRDS(nm, file = "old_nullmod_lmm.rds", version = 2)

# Logistic
nm <- fitNullModel(dat, outcome="outcome_binary", covars="b", verbose=FALSE, family = "binomial")
saveRDS(nm, file = "old_nullmod_glm.rds", version = 2)

# GLMM
nm <- fitNullModel(dat, outcome="outcome_binary", covars="b", verbose=FALSE, family = "binomial", cov.mat = cov_mat)
saveRDS(nm, file = "old_nullmod_glmm.rds", version = 2)

# WLS
nm <- fitNullModel(dat, outcome = "outcome_linear", covars = "b", group.var="group", verbose=FALSE)
saveRDS(nm, file = "old_nullmod_wls.rds", version = 2)
