\name{pcrelateMakeGRM}
\alias{pcrelateMakeGRM}
\title{Creates a Genetic Relationship Matrix (GRM) of Pairwise Kinship Coefficient Estimates from PC-Relate Output}
\description{\code{pcrelateMakeGRM} is used to create a genetic relationship matrix (GRM) of pairwise kinship coefficient estimates from the output of \code{pcrelate}.}
\usage{
pcrelateMakeGRM(pcrelobj, sample.include = NULL, thresh = NULL, scaleKin = 2, verbose = TRUE)
}
\arguments{
  \item{pcrelobj}{The object containing the output from \code{pcrelate}.  This could be a list of class \code{pcrelate} or an object of class \code{gds.class} read into \code{R} using the function \code{\link[gdsfmt]{openfn.gds}} from the \code{\link[gdsfmt:gdsfmt-package]{gdsfmt}} package.}
  \item{sample.include}{A vector of IDs for samples to be included in the GRM. The default is NULL, which includes all samples in \code{pcrelobj}.}
  \item{thresh}{Kinship threshold for clustering samples to make the output matrix sparse block-diagonal. When \code{NULL}, no clustering is done. See 'Details'.}
  \item{scaleKin}{Specifies a numeric constant to scale each estimated kinship coefficient by in the GRM.  The default value is 2.}
  \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is TRUE.}
}
\details{
  This function provides a quick and easy way to construct a genetic relationship matrix (GRM) from the output of \code{pcrelate}.
  
\code{thresh} sets a threhsold for clustering samples such that any pair with an estimated kinship value greater than \code{thresh} is in the same cluster. All pairwise estimates within a cluster are kept, even if they are below \code{thresh}. All pairwise estimates between clusters are set to 0, creating a sparse, block-diagonal matrix. When \code{thresh} is \code{NULL}, no clustering is done and all samples are returned in one block.
}
\value{
An object of class '\code{Matrix}' with pairwise kinship coefficients.
}
\references{
    Conomos M.P., Reiner A.P., Weir B.S., & Thornton T.A. (2016). Model-free Estimation of Recent Genetic Relatedness. American Journal of Human Genetics, 98(1), 127-148.
}
\author{Matthew P. Conomos}
\seealso{
	\code{\link{pcrelate}} for the function that performs PC-Relate.
}
% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
