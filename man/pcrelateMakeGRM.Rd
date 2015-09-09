\name{pcrelateMakeGRM}
\alias{pcrelateMakeGRM}
\title{Creates a Genetic Relationship Matrix (GRM) of Pairwise Kinship Coefficient Estimates from PC-Relate GDS Output}
\description{\code{pcrelateMakeGRM} is used to create a genetic relationship matrix (GRM) of pairwise kinship coefficient and estimates from a GDS file output by \code{pcrelate}.}
\usage{
pcrelateMakeGRM(pcrelObj, scan.include = NULL, scaleKin = 2)
}
\arguments{
  \item{pcrelObj}{The GDS object containing the output from \code{pcrelate}.  This can be created in \code{R} using the function \code{\link[gdsfmt]{openfn.gds}} from the \code{\link[gdsfmt:gdsfmt-package]{gdsfmt}} package.}
  \item{scan.include}{A vector of IDs for samples to be included in the GRM. The default is NULL, which includes all samples in \code{pcrelObj}.}
  \item{scaleKin}{Specifies a numeric constant to scale each estimated kinship coefficient by in the GRM.  The default value is 2.}
}
\details{
	This function is intended to be a helpful utility for those unfamiliar with GDS.  It provides a quick and easy way to construct a genetic relationship matrix (GRM) from the GDS file output by \code{pcrelate}.
}
\references{
    Conomos M.P., Reiner, A.P., Weir, B.S., & Thornton T.A. Model-Free Estimation of Recent Genetic Relatedness. (In Review).
}
\author{Matthew P. Conomos}
\seealso{
	\code{\link{pcrelate}} for the function that performs PC-Relate and creates the GDS output file.
	\code{\link{pcrelateReadKinship}} for the function that creates a table of pairwise kinship coefficient and IBD sharing probabilities from the same GDS output file.
	\code{\link{pcrelateReadInbreed}} for the function that creates a table of inbreeding coefficient estimates from the same GDS output file.
}
% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
