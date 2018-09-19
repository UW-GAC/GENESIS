\name{print.pcair}
\alias{print.pcair}
\alias{summary.pcair}
\alias{print.summary.pcair}
\title{PC-AiR: Principal Components Analysis in Related Samples}
\description{Print methods for pcair}
\usage{
\method{print}{pcair}(x, ...)
\method{summary}{pcair}(object, ...)
\method{print}{summary.pcair}(x, ...)
}
\arguments{
  \item{object}{An object of class '\code{pcair}', i.e. output from the \code{pcair} function. }
  \item{x}{An object of class '\code{pcair}', i.e. output from the \code{pcair} function. }
  \item{\dots}{Further arguments passed to or from other methods. }
}
\author{Matthew P. Conomos}
\seealso{
	\code{\link{pcair}} for obtaining principal components that capture population structure in the presence of relatedness.
}
\examples{
# file path to GDS file
gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
# read in GDS data
HapMap_geno <- gdsfmt::openfn.gds(gdsfile)
# load saved matrix of KING-robust estimates
data("HapMap_ASW_MXL_KINGmat")
# run PC-AiR
mypcair <- pcair(HapMap_geno, kinobj = HapMap_ASW_MXL_KINGmat, 
                divobj = HapMap_ASW_MXL_KINGmat)
print(mypcair)
summary(mypcair)
gdsfmt::closefn.gds(HapMap_geno)
}
\keyword{robust}
\keyword{multivariate}
\keyword{ancestry}

