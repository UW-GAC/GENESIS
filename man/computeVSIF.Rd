\name{computeVSIF}
\alias{computeVSIF}
\title{
Computes variant-specific inflation factors
}
\description{
Computes variant-specific inflation factors resulting from differences 
in variances and allele frequencies across groups pooled togehter in analysis.
}
\usage{
computeVSIF(eafs, ns, sigma_sqs)
}
\arguments{
  \item{eafs}{A named vector of a matrix of effect allele frequencies across
  groups. Vector/column names are group names.}
  \item{ns}{A named vector of group sample sizes}
  \item {sigma_sqs}{A named vector of residual variances across groups}
}
\details{
The function computes the expected inflation/deflation for each specific variant
caused by differences in allele frequencies in combination with differences
in residual variances across groups that are aggregated together (e.g. 
individuals with different genetic ancestry patterns). The inflation/deflation 
is especially expected if a homogeneous variance model is used. 
If a null model object was fit using the GENESIS R package (this package), 
see function \code{\link{computeVSIFnullmod}} to extract sample sizes and
residual variances based on the null model. The null model should be fit under 
a homogeneous variance model.
}
\value{
  \item{SE_true} Large sample test statistic variances accounting for 
  differences in residual variances.
  \item{SE_naive} Large sample test statistic variances (wrongly) assuming 
  that all residual variances are the same across groups.
  \item{Inflation_factor} Variant-specific inflation factors. Values higher
  than 1 suggest inflation (too significant p-value), values lower than 1 
  suggest deflation (too high p-value).
}
\references{
Sofer, T., Zheng, X., Laurie, C. A., Gogarten, S. M., Brody, J. A., 
Conomos, M. P., ... & Rice, K. M. (2020). Population Stratification at the 
Phenotypic Variance level and Implication for the Analysis of 
Whole Genome Sequencing Data from Multiple Studies. BioRxiv.
}
\author{
Tamar Sofer, Kenneth M Rice
}
\examples{
  ns <- c(2000, 5000, 100)
  sigma_sqs <- c(1, 1, 2)
  eaf_vec <- c(0.1, 0.2, 0.5)
  names(eaf_vec) <- names(ns) <- names(sigma_sqs) <- c("g1", "g2", "g3") 
  res <- computeVSIF(eafs = eaf_vec, ns, sigma_sqs)
  
  eafs_mat <- matrix(c(0.1, 0.2, 0.5, 0.1, 0.01, 0.5), nrow = 2, byrow = TRUE)
  colnames(eafs_mat) <- names(sigma_sqs)
  res <- computeVSIF(eafs = eafs_mat, ns, sigma_sqs)
}

