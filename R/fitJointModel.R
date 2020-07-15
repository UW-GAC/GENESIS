# Inputs:
# * nullmod: null model object without variants
# * G: genotype matrix with variants to estimate jointly

# Outputs:
# * list with elements: pve, fixef, covar
jointScoreTest <- function(nullmod, G) {  # # Check rownames/colnames match.

  # Check that all samples in null model are in the genotype matrix.
  missing_ids <- setdiff(nullmod$sample.id, rownames(G))
  if (length(missing_ids) > 0) {
    stop("missing samples in genotype matrix!")
  }

  # Reorder samples.
  idx <- match(nullmod$sample.id, rownames(G))
  G <- G[idx, , drop = FALSE]

  # Genotype adjusted for covariates and random effects.
  Gtilde <- calcGtilde(nullmod, G)

  # Score statistic.
  GY <- crossprod(Gtilde, nullmod$Ytilde)

  # GPG.
  GG <- crossprod(Gtilde)

  # Joint test statistic. Note this is the equivalent of the stat being squared, by convention.
  Stat.joint <- crossprod(GY, solve(GG, GY))
  pval.joint <- pchisq(Stat.joint, df = ncol(G))

  # Calculate proportion of variance explained.
  pve <- as.vector(Stat.joint / nullmod$RSS0)

  # Covariance matrix.
  betaCov <- solve(GG)
  rownames(betaCov) <- colnames(betaCov) <- colnames(G)

  beta <- as.vector(solve(GG, GY))
  se <- sqrt(diag(betaCov))

  # Calculate fixed effects data frame.
  fixef <- data.frame(
    Est = beta,
    SE = se,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      Stat = beta / se,
      pval = pchisq(Stat^2, lower.tail = F, df = 1)
    )
  rownames(fixef) <- colnames(G)

  res <- list()

  res$Stat.joint <- as.numeric(Stat.joint)
  res$pval.joint <- as.numeric(pval.joint)
  res$pve <- pve
  res$fixef <- fixef
  res$betaCov <- betaCov

  return(res)

}
