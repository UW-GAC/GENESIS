# Inputs:
# * nullmod: null model object without variants
# * G: genotype matrix with variants to estimate jointly

# Outputs:
# * list with elements: pve, fixef?

.fitJointModel <- function(nullmod, G) {  # # Check rownames/colnames match.

  res <- list()

  res$pve <- numeric()
  res$fixef <- data.frame(stringsAsFactors = FALSE)
  res$covar <- matrix()

  return(res)

}
