\name{pcrelateReadKinship}
\alias{pcrelateReadKinship}
\title{Create a Table of Pairwise Kinship Coefficient and IBD Sharing Probability Estimates from PC-Relate Output}
\description{\code{pcrelateReadKinship} is used to create a table of pairwise kinship coefficient and IBD sharing probability (k0, k1, k2) estimates from the output of \code{pcrelate}.}
\usage{
pcrelateReadKinship(pcrelObj, scan.include = NULL, ibd.probs = TRUE, 
            kin.thresh = NULL)
}
\arguments{
  \item{pcrelObj}{The object containing the output from \code{pcrelate}.  This could be a list of class \code{pcrelate} or an object of class \code{gds.class} read into \code{R} using the function \code{\link[gdsfmt]{openfn.gds}} from the \code{\link[gdsfmt:gdsfmt-package]{gdsfmt}} package.}
  \item{scan.include}{A vector of IDs for samples to be included in the table. The default is NULL, which includes all samples in \code{pcrelObj}.}
  \item{ibd.probs}{Logical indicator of whether or not the output in \code{pcrelObj} has estimates of IBD sharing probabilities.}
  \item{kin.thresh}{Specifies a minimum value of the estimated kinship coefficient to include in the table; i.e. only pairs with an estimated kinship coefficient greater than \code{kin.thresh} will be included in the table.  The default is NULL, which includes all pairs.}
}
\details{
	This function provides an easy way to make a table of pairwise relatedness estimates.
}
\references{
    Conomos M.P., Reiner A.P., Weir B.S., & Thornton T.A. (2016). Model-free Estimation of Recent Genetic Relatedness. American Journal of Human Genetics, 98(1), 127-148.
}
\author{Matthew P. Conomos}
\seealso{
	\code{\link{pcrelate}} for the function that performs PC-Relate.
	\code{\link{pcrelateReadInbreed}} for the function that creates a table of inbreeding coefficient estimates from the same PC-Relate output file.
	\code{\link{pcrelateMakeGRM}} for the function that creates a genetic relationship matrix (GRM) of pairwise kinship coefficient estimates from the same PC-Relate output file.
}
% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
