# Inputs:
# * nullmod: null model object without variants
# * G: genotype matrix with variants to estimate jointly

# Outputs:
# * list with elements: pve, fixef, covar
jointScoreTest <- function(nullmod, G) {  # # Check rownames/colnames match.

  # Check that all samples in null model are in the genotype matrix.
  missing_ids <- setdiff(rownames(nullmod$model.matrix), rownames(G))
  if (length(missing_ids) > 0) {
    stop("missing samples in genotype matrix!")
  }

  # Check for missingness in genotype data.
  # For now, missingness is not allowed.
  # Eventually, update to handle missingness the same way as assocTestSingle.
  if (any(is.na(G))) {
    stop("genotype matrix cannot contain missing data!")
  }

  # Reorder samples.
  idx <- match(rownames(nullmod$model.matrix), rownames(G))
  G <- G[idx, , drop = FALSE]

  # Genotype adjusted for covariates and random effects.
  Gtilde <- calcGtilde(nullmod, G)

  # Score statistic.
  GY <- crossprod(Gtilde, nullmod$Ytilde)

  # GPG.
  GG <- crossprod(Gtilde)

  # Covariance matrix for estimates.
  betaCov <- solve(GG)
  rownames(betaCov) <- colnames(betaCov) <- colnames(G)

  # Fixed effect estimates for the variants.
  beta <- crossprod(betaCov, GY)
  se <- sqrt(diag(betaCov))

  # Joint test statistic. Note this is the equivalent of the stat being squared,
  # by convention.
  Stat.joint <- as.numeric(crossprod(GY, beta))
  pval.joint <- pchisq(Stat.joint, lower.tail = FALSE, df = ncol(G))

  # Percentage of variance explained jointly by these variants.
  pve.joint <- as.numeric(Stat.joint / nullmod$RSS0)

  # Create fixed effects data frame.
  fixef <- data.frame(
    Est = as.vector(beta),
    Est.SE = se,
    stringsAsFactors = FALSE
  )
  fixef$Stat <- fixef$Est / fixef$Est.SE
  fixef$pval <- pchisq(fixef$Stat^2, lower.tail = FALSE, df = 1)
  fixef$PVE <- fixef$Stat^2 / nullmod$RSS0

  rownames(fixef) <- colnames(G)

  res <- list()

  res$Joint.Stat <- as.numeric(Stat.joint)
  res$Joint.pval <- as.numeric(pval.joint)
  res$Joint.PVE <- pve.joint
  res$fixef <- fixef
  res$betaCov <- betaCov

  return(res)

}
