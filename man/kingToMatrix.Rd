\name{kingToMatrix}
\alias{kingToMatrix}
\alias{kingToMatrix,character-method}
\alias{kingToMatrix,snpgdsIBDClass-method}
\title{Convert KING text output to an R Matrix}
\description{\code{kingToMatrix} is used to extract the pairwise kinship coefficient estimates from the output text files of KING --ibdseg, KING --kinship, or KING --related and put them into an R object of class \code{Matrix}. One use of this matrix is that it can be read by the functions \code{\link{pcair}} and \code{\link{pcairPartition}}.}
\usage{
\S4method{kingToMatrix}{character}(king, estimator = c("PropIBD", "Kinship"), sample.include = NULL, thresh = NULL, verbose = TRUE)
\S4method{kingToMatrix}{snpgdsIBDClass}(king, sample.include = NULL, thresh = 2^(-11/2), verbose = TRUE)
}
\arguments{
  \item{king}{Output from KING, either a \code{snpgdsIBDClass} object from \code{\link{snpgdsIBDKING}} or a character vector of one or more file names output from the command-line version of KING; see 'Details'.}
  \item{estimator}{Which estimates to read in when using command-line KING output; must be either "PropIBD" or "Kinship"; see 'Details'.}
  \item{sample.include}{An optional vector of sample.id indicating all samples that should be included in the output matrix; see 'Details' for usage.}
  \item{thresh}{Kinship threshold for clustering samples to make the output matrix sparse block-diagonal. When \code{NULL}, no clustering is done. See 'Details'.}
  \item{verbose}{A logical indicating whether or not to print status updates to the console; the default is TRUE.}
}
\details{\code{king} can be a vector of multiple file names if your KING output is stored in multiple files; e.g. KING --kinship run with family IDs returns a .kin and a .kin0 file for pairs within and not within the same family, respectively.

When reading command-line KING output, the \code{estimator} argument is required to specify which estimates to read in. When reading KING --ibdseg output, only "PropIBD" will be available; when reading KING --kinship output, only "Kinship" will be available; when reading KING --related output, both "PropIBD" and "Kinship" will be available - use this argument to select which to read. See the KING documentation for details on each estimator. 

\code{sample.include} has two primary functions: 1) It can be used to subset the KING output. 2) \code{sample.include} can include sample.id not in \code{king}; this ensures that all samples will be in the output matrix when reading KING --ibdseg output, which likely does not contain all pairs. When left \code{NULL}, the function will determine the list of samples from what is observed in \code{king}. It is recommended to use \code{sample.include} to ensure all of your samples are included in the output matrix.

\code{thresh} sets a threhsold for clustering samples such that any pair with an estimated kinship value greater than \code{thresh} is in the same cluster. All pairwise estimates within a cluster are kept, even if they are below \code{thresh}. All pairwise estimates between clusters are set to 0, creating a sparse, block-diagonal matrix. When \code{thresh} is \code{NULL}, no clustering is done and all samples are returned in one block. This feature is useful when converting KING --ibdseg or KING --robust estimates to be used as a kinship matrix, if you have a lower threshold that you consider 'related'. This feature should not be used when converting KING --robust estimates to be used as \code{divobj} in \code{\link{pcair}} or \code{\link{pcairPartition}}, as PC-AiR requires the negative estimates to identify ancestrally divergent pairs.}
\value{An object of class '\code{Matrix}' with pairwise kinship coefficients by KING --ibdseg or KING --robust for each pair of individuals in the sample.  The estimates are on both the upper and lower triangle of the matrix, and the diagonal is arbitrarily set to 0.5. sample.id are set as the column and row names of the matrix.}
\references{
    Conomos M.P., Miller M., & Thornton T. (2015). Robust Inference of Population Structure for Ancestry Prediction and Correction of Stratification in the Presence of Relatedness. Genetic Epidemiology, 39(4), 276-293.
    
    Manichaikul, A., Mychaleckyj, J.C., Rich, S.S., Daly, K., Sale, M., & Chen, W.M. (2010). Robust relationship inference in genome-wide association studies. Bioinformatics, 26(22), 2867-2873.
}
\author{Matthew P. Conomos}
\seealso{
	\code{\link{pcair}} and \code{\link{pcairPartition}} for functions that use the output matrix.
}
\examples{
# KING --kinship
file.king <- c(system.file("extdata", "MXL_ASW.kin0", package="GENESIS"),
               system.file("extdata", "MXL_ASW.kin", package="GENESIS"))
KINGmat <- kingToMatrix(file.king, estimator="Kinship")

# KING --ibdseg
file.king <- system.file("extdata", "HapMap.seg", package="GENESIS")
KINGmat <- kingToMatrix(file.king, estimator="PropIBD")

# SNPRelate
library(SNPRelate)
gds <- snpgdsOpen(system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS"))
king <- snpgdsIBDKING(gds)
KINGmat <- kingToMatrix(king)
snpgdsClose(gds)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
