# Inputs:
# * null.model: null model object without variants
# * G: genotype matrix with variants to estimate jointly

# Outputs:
# * list with elements: pve, fixef, covar
jointScoreTest <- function(null.model, G) {  # # Check rownames/colnames match.

  # Convert old null model format if necessary.
  null.model <- .updateNullModelFormat(null.model)

  # check that null model has required elements
  if (isNullModelSmall(null.model)) {
    stop("small null model cannot be used for a joint score test")
  }

  # Check that all samples in null model are in the genotype matrix.
  missing_ids <- setdiff(rownames(null.model$model.matrix), rownames(G))
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
  idx <- match(rownames(null.model$model.matrix), rownames(G))
  G <- G[idx, , drop = FALSE]

  # Genotype adjusted for covariates and random effects.
  Gtilde <- calcGtilde(null.model, G)

  # Score statistic.
  GY <- crossprod(Gtilde, null.model$fit$resid.cholesky)

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
  pve.joint <- as.numeric(Stat.joint / null.model$RSS0)

  # Create fixed effects data frame.
  fixef <- data.frame(
    Est = as.vector(beta),
    SE = se,
    stringsAsFactors = FALSE
  )
  fixef$Stat <- fixef$Est / fixef$SE
  fixef$pval <- pchisq(fixef$Stat^2, lower.tail = FALSE, df = 1)
  fixef$PVE <- fixef$Stat^2 / null.model$RSS0

  rownames(fixef) <- colnames(G)

  res <- list()

  res$Joint.Stat <- as.numeric(Stat.joint)
  res$Joint.pval <- as.numeric(pval.joint)
  res$Joint.PVE <- pve.joint
  res$fixef <- fixef
  res$betaCov <- betaCov

  return(res)

}
