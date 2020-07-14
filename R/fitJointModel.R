# Inputs:
# * nullmod: null model object without variants
# * G: genotype matrix with variants to estimate jointly

# Outputs:
# * list with elements: pve, fixef?
fitJointModel <- function(nullmod, G) {  # # Check rownames/colnames match.

  # Genotype adjusted for covariates and random effects.
  Gtilde <- calcGtilde(nullmod, G)

  # Score statistic.
  GY <- crossprod(Gtilde, nullmod$Ytilde)

  # GPG.
  GG <- crossprod(Gtilde)

  # Calculate proportion of variance explained.
  pve <- as.vector(crossprod(GY, solve(GG, GY)) / nullmod$RSS0)

  # Covariance matrix.
  covar <- solve(GG)
  rownames(covar) <- colnames(covar) <- colnames(G)

  beta <- as.vector(solve(GG, GY))
  se <- sqrt(diag(covar))

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

  res$pve <- pve
  res$fixef <- fixef
  res$covar <- covar

  return(res)

}
