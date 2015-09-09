\name{pcrelateReadKinship}
\alias{pcrelateReadKinship}
\title{Create a Table of Pairwise Kinship Coefficient and IBD Sharing Probability Estimates from PC-Relate GDS Output}
\description{\code{pcrelateReadKinship} is used to create a table of pairwise kinship coefficient and IBD sharing probability (k0, k1, k2) estimates from a GDS file output by \code{pcrelate}.}
\usage{
pcrelateReadKinship(pcrelObj, scan.include = NULL, ibd.probs = TRUE, 
            kin.thresh = NULL)
}
\arguments{
  \item{pcrelObj}{The GDS object containing the output from \code{pcrelate}.  This can be created in \code{R} using the function \code{\link[gdsfmt]{openfn.gds}} from the \code{\link[gdsfmt:gdsfmt-package]{gdsfmt}} package.}
  \item{scan.include}{A vector of IDs for samples to be included in the table. The default is NULL, which includes all samples in \code{pcrelObj}.}
  \item{ibd.probs}{Logical indicator of whether or not the output in \code{pcrelObj} has estimates of IBD sharing probabilities.}
  \item{kin.thresh}{Specifies a minimum value of the estimated kinship coefficient to include in the table; i.e. only pairs with an estimated kinship coefficient greater than \code{kin.thresh} will be included in the table.  The default is NULL, which includes all pairs.}
}
\details{
	This function is intended to be a helpful utility for those unfamiliar with GDS.
}
\references{
    Conomos M.P., Reiner, A.P., Weir, B.S., & Thornton T.A. Model-Free Estimation of Recent Genetic Relatedness. (In Review).
}
\author{Matthew P. Conomos}
\seealso{
	\code{\link{pcrelate}} for the function that performs PC-Relate and creates the GDS output file.
	\code{\link{pcrelateReadInbreed}} for the function that creates a table of inbreeding coefficient estimates from the same GDS output file.
	\code{\link{pcrelateMakeGRM}} for the function that creates a genetic relationship matrix (GRM) of pairwise kinship coefficient estimates from the same GDS output file.
}
% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
