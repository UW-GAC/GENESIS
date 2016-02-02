\name{pcrelateReadInbreed}
\alias{pcrelateReadInbreed}
\title{Create a Table of Inbreeding Coefficient Estimates from PC-Relate Output}
\description{\code{pcrelateReadInbreed} is used to create a table of inbreeding coefficient estimates from the output of \code{pcrelate}.}
\usage{
pcrelateReadInbreed(pcrelObj, scan.include = NULL, f.thresh = NULL)
}
\arguments{
  \item{pcrelObj}{The object containing the output from \code{pcrelate}.  This could be a list of class \code{pcrelate} or an object of class \code{gds.class} read into \code{R} using the function \code{\link[gdsfmt]{openfn.gds}} from the \code{\link[gdsfmt:gdsfmt-package]{gdsfmt}} package.}
  \item{scan.include}{A vector of IDs for samples to be included in the table. The default is NULL, which includes all samples in \code{pcrelObj}.}  
  \item{f.thresh}{Specifies a minimum value of the estimated inbreeding coefficient to include in the table; i.e. only individuals with an estimated inbreeding coefficient greater than \code{f.thresh} will be included in the table.  The default is NULL, which includes all individuals.}
}
\details{
	This function provides an easy way to make a table of estimated inbreeding coefficients.
}
\references{
    Conomos M.P., Reiner A.P., Weir B.S., & Thornton T.A. (2016). Model-free Estimation of Recent Genetic Relatedness. American Journal of Human Genetics, 98(1), 127-148.
}
\author{Matthew P. Conomos}
\seealso{
	\code{\link{pcrelate}} for the function that performs PC-Relate.
	\code{\link{pcrelateReadKinship}} for the function that creates a table of pairwise kinship coefficient and IBD sharing probabilities from the same PC-Relate output file.
	\code{\link{pcrelateMakeGRM}} for the function that creates a genetic relationship matrix (GRM) of pairwise kinship coefficient estimates from the same PC-Relate output file.
}
% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
