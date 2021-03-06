\name{pcairPartition}
\alias{pcairPartition}
\title{Partition a sample into an ancestry representative 'unrelated subset' and a 'related subset'}
\description{\code{pcairPartition} is used to partition a sample from a genetic study into an ancestry representative 'unrelated subset' and a 'related subset'.  The 'unrelated subset' contains individuals who are all mutually unrelated to each other and representative of the ancestries of all individuals in the sample, and the 'related subset' contains individuals who are related to someone in the 'unrealted subset'.}
\usage{
pcairPartition(kinobj, divobj = NULL,
               kin.thresh = 2^(-11/2), div.thresh = -2^(-11/2),
               unrel.set = NULL, sample.include = NULL, verbose = TRUE)
}
\arguments{
  \item{kinobj}{A symmetric matrix of pairwise kinship coefficients for every pair of individuals in the sample: upper and lower triangles must both be filled; diagonals should be self-kinship or set to a non-missing constant value. This matrix is used for partitioning the sample into the 'unrelated' and 'related' subsets.  See 'Details' for how this interacts with \code{kin.thresh} and \code{unrel.set}.  IDs for each individual must be set as the column names of the matrix. This matrix may also be provided as a GDS object; see 'Details'.}
  \item{divobj}{A symmetric matrix of pairwise ancestry divergence measures for every pair of individuals in the sample: upper and lower triangles must both be filled; diagonals should be set to a non-missing constant value. This matrix is used for partitioning the sample into the 'unrelated' and 'related' subsets.  See 'Details' for how this interacts with \code{div.thresh}.  IDs for each individual must be set as the column names of the matrix.This matrix may be identical to \code{kinobj}. This matrix may be \code{NULL} to ignore ancestry divergence. This matrix may also be provided as a GDS object; see 'Details'.}
  \item{kin.thresh}{Threshold value on \code{kinobj} used for declaring each pair of individuals as related or unrelated.  The default value is 2^(-11/2) ~ 0.022, corresponding to 4th degree relatives.  See 'Details' for how this interacts with \code{kinobj}.}
  \item{div.thresh}{Threshold value on \code{divobj} used for deciding if each pair of individuals is ancestrally divergent.  The default value is -2^(-11/2) ~ -0.022.  See 'Details' for how this interacts with \code{divobj}.}
  \item{unrel.set}{An optional vector of IDs for identifying individuals that are forced into the unrelated subset.  See 'Details' for how this interacts with \code{kinobj}.}
  \item{sample.include}{An optional vector of IDs for selecting samples to consider for either set.}
  \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is TRUE.}
}
\details{
We recommend using software that accounts for population structure to estimate pairwise kinship coefficients to be used in \code{kinobj}.  Any pair of individuals with a pairwise kinship greater than \code{kin.thresh} will be declared 'related.'  Kinship coefficient estimates from the KING-robust software are typically used as measures of ancestry divergence in \code{divobj}.  Any pair of individuals with a pairwise divergence measure less than \code{div.thresh} will be declared ancestrally 'divergent'.  Typically, \code{kin.thresh} and \code{div.thresh} are set to be the amount of error around 0 expected in the estimate for a pair of truly unrelated individuals.  If \code{unrel.set = NULL}, the PC-AiR algorithm is used to find an 'optimal' partition (see 'References' for a paper describing the algorithm).  If \code{unrel.set} and \code{kinobj} are both specified, then all individuals with IDs in \code{unrel.set} are forced in the 'unrelated subset' and the PC-AiR algorithm is used to partition the rest of the sample; this is especially useful for including reference samples of known ancestry in the 'unrelated subset'.

For large sample sizes, storing both \code{kinobj} and \code{divobj} in memory may be prohibitive. Both matrices may be stored in GDS files and provided as \code{gds.class} objects. \code{\link{mat2gds}} saves matrices in GDS format. Alternatively, \code{kinobj} (but not \code{divobj}) can be represented as a sparse \code{\link{Matrix}} object; see \code{\link{kingToMatrix}} and \code{\link{pcrelateToMatrix}}.

Matrix objects from the \pkg{\link{Matrix}} package are also supported.
}
\value{A list including:
	\item{rels}{A vector of IDs for individuals in the 'related subset'.}
	\item{unrels}{A vector of IDs for individuals in the 'unrelated subset'.}
}
\references{
	Conomos M.P., Miller M., & Thornton T. (2015). Robust Inference of Population Structure for Ancestry Prediction and Correction of Stratification in the Presence of Relatedness. Genetic Epidemiology, 39(4), 276-293.
    
    Manichaikul, A., Mychaleckyj, J.C., Rich, S.S., Daly, K., Sale, M., & Chen, W.M. (2010). Robust relationship inference in genome-wide association studies. Bioinformatics, 26(22), 2867-2873.
}
\author{Matthew P. Conomos}
\note{\code{pcairPartition} is called internally in the function \code{pcair} but may also be used on its own to partition the sample into an ancestry representative 'unrelated' subset and a 'related' subset without performing PCA.}

\seealso{
	\code{\link{pcair}} which uses this function for finding principal components in the presence of related individuals.
	\code{\link{kingToMatrix}} for creating a matrix of kinship coefficent estimates or pairwise ancestry divergence measures from KING output text files that can be used as \code{kinobj} or \code{divobj}.
	\code{\link{kin2gds}} and \code{\link{mat2gds}} for saving kinship matrices to GDS.
}
\examples{
# load saved matrix of KING-robust estimates
data("HapMap_ASW_MXL_KINGmat")
# partition the sample
part <- pcairPartition(kinobj = HapMap_ASW_MXL_KINGmat, 
                       divobj = HapMap_ASW_MXL_KINGmat)
}

