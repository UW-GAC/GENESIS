\name{pcrelateReadInbreed}
\alias{pcrelateReadInbreed}
\title{Create a Table of Inbreeding Coefficient Estimates from PC-Relate GDS Output}
\description{\code{pcrelateReadInbreed} is used to create a table of inbreeding coefficient estimates from a GDS file output by \code{pcrelate}.}
\usage{
pcrelateReadInbreed(pcrelObj, scan.include = NULL, f.thresh = NULL)
}
\arguments{
  \item{pcrelObj}{The GDS object containing the output from \code{pcrelate}.  This can be created in \code{R} using the function \code{\link[gdsfmt]{openfn.gds}} from the \code{\link[gdsfmt:gdsfmt-package]{gdsfmt}} package.}
  \item{scan.include}{A vector of IDs for samples to be included in the table. The default is NULL, which includes all samples in \code{pcrelObj}.}  
  \item{f.thresh}{Specifies a minimum value of the estimated inbreeding coefficient to include in the table; i.e. only individuals with an estimated inbreeding coefficient greater than \code{f.thresh} will be included in the table.  The default is NULL, which includes all individuals.}
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
	\code{\link{pcrelateReadKinship}} for the function that creates a table of pairwise kinship coefficient and IBD sharing probabilities from the same GDS output file.
	\code{\link{pcrelateMakeGRM}} for the function that creates a genetic relationship matrix (GRM) of pairwise kinship coefficient estimates from the same GDS output file.
}
% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
