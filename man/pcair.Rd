\name{pcair}
\alias{pcair}
\alias{pcair-methods}
\alias{pcair,gds.class-method}
\alias{pcair,SNPGDSFileClass-method}
\alias{pcair,GdsGenotypeReader-method}
\alias{pcair,GenotypeData-method}
\alias{pcair,SeqVarGDSClass-method}
\title{PC-AiR: Principal Components Analysis in Related Samples}
\description{\code{pcair} is used to perform a Principal Components Analysis using genome-wide SNP data for the detection of population structure in a sample.  Unlike a standard PCA, PC-AiR accounts for sample relatedness (known or cryptic) to provide accurate ancestry inference that is not confounded by family structure.}
\usage{
\S4method{pcair}{gds.class}(gdsobj, kinobj, divobj,
                            kin.thresh = 2^(-11/2), div.thresh = -2^(-11/2),
                            unrel.set = NULL, 
                            sample.include = NULL, snp.include = NULL,
                            num.thread = 1L, verbose = TRUE, ...)
\S4method{pcair}{SNPGDSFileClass}(gdsobj, ...)
\S4method{pcair}{GdsGenotypeReader}(gdsobj, ...)
\S4method{pcair}{GenotypeData}(gdsobj, ...)
\S4method{pcair}{SeqVarGDSClass}(gdsobj, ...)
}
\arguments{
  \item{gdsobj}{An object providing a connection to a GDS file.}
  \item{kinobj}{A symmetric matrix of pairwise kinship coefficients for every pair of individuals in the sample: upper and lower triangles must both be filled; diagonals should be self-kinship or set to a non-missing constant value. This matrix is used for partitioning the sample into the 'unrelated' and 'related' subsets.  See 'Details' for how this interacts with \code{kin.thresh} and \code{unrel.set}.  IDs for each individual must be set as the column names of the matrix. This matrix may also be provided as a GDS object; see 'Details'.}
  \item{divobj}{A symmetric matrix of pairwise ancestry divergence measures for every pair of individuals in the sample: upper and lower triangles must both be filled; diagonals should be set to a non-missing constant value. This matrix is used for partitioning the sample into the 'unrelated' and 'related' subsets.  See 'Details' for how this interacts with \code{div.thresh}.  IDs for each individual must be set as the column names of the matrix. This matrix may be identical to \code{kinobj}. This matrix may also be provided as a GDS object; see 'Details'.}
  \item{kin.thresh}{Threshold value on \code{kinobj} used for declaring each pair of individuals as related or unrelated.  The default value is 2^(-11/2) ~ 0.02, corresponding to 4th degree relatives.  See 'Details' for how this interacts with \code{kinobj}.}
  \item{div.thresh}{Threshold value on \code{divobj} used for deciding if each pair of individuals is ancestrally divergent.  The default value is -2^(-11/2) ~ -0.022.  See 'Details' for how this interacts with \code{divobj}.}
  \item{unrel.set}{An optional vector of IDs for identifying individuals that are forced into the unrelated subset.  See 'Details' for how this interacts with \code{kinobj}.}
  \item{sample.include}{An optional vector of IDs for selecting samples to consider for either set.}
  \item{snp.include}{An optional vector of snp or variant IDs to use in the analysis.}
  \item{num.thread}{Number of threads; passed to \code{\link{snpgdsPCA}}.}
  \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is TRUE.}
  \item{\dots}{Additional arguments to pass to \code{\link{snpgdsPCA}}.}
}
\details{The basic premise of PC-AiR is to partition the entire sample of individuals into an ancestry representative 'unrelated subset' and a 'related set', perform standard PCA on the 'unrelated subset', and predict PC values for the 'related subset'.
	
	We recommend using software that accounts for population structure to estimate pairwise kinship coefficients to be used in \code{kinobj}.  Any pair of individuals with a pairwise kinship greater than \code{kin.thresh} will be declared 'related.'  Kinship coefficient estimates from the KING-robust software are typically used as measures of ancestry divergence in \code{divobj}.  Any pair of individuals with a pairwise divergence measure less than \code{div.thresh} will be declared ancestrally 'divergent'.  Typically, \code{kin.thresh} and \code{div.thresh} are set to be the amount of error around 0 expected in the estimate for a pair of truly unrelated individuals.
	
	\code{kinobj} and \code{divobj} may be identical.
		
	There are multiple ways to partition the sample into an ancestry representative 'unrelated subset' and a 'related subset'.  If \code{unrel.set = NULL}, then the PC-AiR algorithm is used to find an 'optimal' partition (see 'References' for a paper describing the algorithm).  If \code{unrel.set} is specified, then all individuals with IDs in \code{unrel.set} are forced in the 'unrelated subset' and the PC-AiR algorithm is used to partition the rest of the sample; this is especially useful for including reference samples of known ancestry in the 'unrelated subset'.
}
\value{An object of class '\code{pcair}'.  A list including:
	\item{vectors}{A matrix of principal components; each column is a principal component. Sample IDs are provided as rownames. The number of PCs returned can be adjusted by supplying the \code{eigen.cnt} argument, which is passed to \code{\link{snpgdsPCA}}.}
	\item{values}{A vector of eigenvalues matching the principal components.  These values are determined from the standard PCA run on the 'unrelated subset'.}
	\item{rels}{A vector of IDs for individuals in the 'related subset'.}
	\item{unrels}{A vector of IDs for individuals in the 'unrelated subset'.}
	\item{kin.thresh}{The threshold value used for declaring each pair of individuals as related or unrelated.}
	\item{div.thresh}{The threshold value used for determining if each pair of individuals is ancestrally divergent.}
	\item{nsamp}{The total number of samples in the analysis.}
	\item{nsnps}{The total number of SNPs used in the analysis.}
	\item{varprop}{The variance proportion for each principal component.}
	\item{call}{The function call passed to \code{pcair}.}
	\item{method}{A character string.  Either "PC-AiR" or "Standard PCA" identifying which method was used for computing principal components. "Standard PCA" is used if no relatives were identified in the sample.}
}
\references{
	Conomos M.P., Miller M., & Thornton T. (2015). Robust Inference of Population Structure for Ancestry Prediction and Correction of Stratification in the Presence of Relatedness. Genetic Epidemiology, 39(4), 276-293.
    
    Manichaikul, A., Mychaleckyj, J.C., Rich, S.S., Daly, K., Sale, M., & Chen, W.M. (2010). Robust relationship inference in genome-wide association studies. Bioinformatics, 26(22), 2867-2873.
}
\author{Matthew P. Conomos}
\seealso{
	\code{\link{pcairPartition}} for a description of the function used by \code{pcair} that can be used to partition the sample into 'unrelated' and 'related' subsets without performing PCA.
	\code{\link{plot.pcair}} for plotting.
	\code{\link{king2mat}} for creating a matrix of pairwise kinship coefficient estimates from KING output text files that can be used for \code{kinobj} or \code{divobj}.
	\code{\link[GWASTools:GWASTools-package]{GWASTools}} for a description of the package containing the following functions: \code{\link{GenotypeData}} for a description of creating a \code{GenotypeData} class object for storing sample and SNP genotype data, \code{\link{MatrixGenotypeReader}} for a description of reading in genotype data stored as a matrix, and \code{\link{GdsGenotypeReader}} for a description of reading in genotype data stored as a GDS file.  Also see \code{\link[SNPRelate]{snpgdsBED2GDS}} in the \code{\link[SNPRelate:SNPRelate-package]{SNPRelate}} package for a description of converting binary PLINK files to GDS.  The generic functions \code{\link{summary}} and \code{\link{print}}.
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
gdsfmt::closefn.gds(HapMap_geno)
}
\keyword{robust}
\keyword{multivariate}
\keyword{ancestry}

