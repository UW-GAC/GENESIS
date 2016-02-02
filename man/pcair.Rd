\name{pcair}
\alias{pcair}
\alias{print.pcair}
\alias{summary.pcair}
\alias{print.summary.pcair}
\title{PC-AiR: Principal Components Analysis in Related Samples}
\description{\code{pcair} is used to perform a Principal Components Analysis using genome-wide SNP data for the detection of population structure in a sample.  Unlike a standard PCA, PC-AiR accounts for sample relatedness (known or cryptic) to provide accurate ancestry inference that is not confounded by family structure.}
\usage{
pcair(genoData, v = 20, kinMat = NULL, kin.thresh = 2^(-11/2), 
            divMat = NULL, div.thresh = -2^(-11/2), unrel.set = NULL, 
            scan.include = NULL, snp.include = NULL, chromosome = NULL,
            snp.block.size = 10000, MAF = 0.01, verbose = TRUE)
\method{print}{pcair}(x, ...)
\method{summary}{pcair}(object, ...)
\method{print}{summary.pcair}(x, ...)
}
\arguments{
  \item{genoData}{An object of class \code{GenotypeData} from the package \code{GWASTools} containing the genotype data for SNPs and samples to be used for the analysis.  This object can easily be created from a matrix of SNP genotype data, PLINK files, or GDS files.}
  \item{v}{The number of principal components to be returned; the default is 20.  If \code{v = NULL}, then all the principal components are returned.}  
  \item{kinMat}{An optional symmetric matrix of pairwise kinship coefficients for every pair of individuals in the sample (the values on the diagonal do not matter, but the upper and lower triangles must both be filled) used for partitioning the sample into the 'unrelated' and 'related' subsets.  See 'Details' for how this interacts with \code{kin.thresh} and \code{unrel.set}.  IDs for each individual must be set as the row and column names of the matrix.}
  \item{kin.thresh}{Threshold value on \code{kinMat} used for declaring each pair of individuals as related or unrelated.  The default value is 2^(-11/2) ~ 0.022.  See 'Details' for how this interacts with \code{kinMat}.}
  \item{divMat}{An optional symmetric matrix of pairwise divergence measures for every pair of individuals in the sample (the values on the diagonal do not matter, but the upper and lower triangles must both be filled) used for partitioning the sample into the 'unrelated' and 'related' subsets.  See 'Details' for how this interacts with \code{div.thresh}.  IDs for each individual must be set as the row and column names of the matrix.}
  \item{div.thresh}{Threshold value on \code{divMat} used for deciding if each pair of individuals is ancestrally divergent.  The default value is -2^(-11/2) ~ -0.022.  See 'Details' for how this interacts with \code{divMat}.}
  \item{unrel.set}{An optional vector of IDs for identifying individuals that are forced into the unrelated subset.  See 'Details' for how this interacts with \code{kinMat}.}
  \item{scan.include}{A vector of IDs for samples to include in the analysis.  If NULL, all samples are included.}
  \item{snp.include}{A vector of SNP IDs to include in the analysis.  If NULL, see \code{chromosome} for further details.}
  \item{chromosome}{A vector of integers specifying which chromosomes to analyze.  This parameter is only considred when \code{snp.include} is NULL; if \code{chromosome} is also NULL, then all SNPs are included.}
  \item{snp.block.size}{The number of SNPs to read-in/analyze at once.  The default value is 10000. }
  \item{MAF}{Minor allele frequency filter; any SNPs with MAF less than this value will be excluded from the analysis; the default value is 0.01.}
  \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is TRUE.}
  \item{object}{An object of class '\code{pcair}', i.e. output from the \code{pcair} function. }
  \item{x}{An object of class '\code{pcair}', i.e. output from the \code{pcair} function. }
  \item{\dots}{Further arguments passed to or from other methods. }
}
\details{The basic premise of PC-AiR is to partition the entire sample of individuals into an ancestry representative 'unrelated subset' and a 'related set', perform standard PCA on the 'unrelated subset', and predict PC values for the 'related subset'.
	
	We recommend using software that accounts for population structure to estimate pairwise kinship coefficients to be used in \code{kinMat}.  Any pair of individuals with a pairwise kinship greater than \code{kin.thresh} will be declared 'related.'  Kinship coefficient estimates from the KING-robust software are used as measures of ancestry divergence in \code{divMat}.  Any pair of individuals with a pairwise divergence measure less than \code{div.thresh} will be declared ancestrally 'divergent'.  Typically, \code{kin.thresh} and \code{div.thresh} are set to be the amount of error around 0 expected in the estimate for a pair of truly unrelated individuals.
	
	If \code{divMat = NULL} and \code{kinMat} is specified, the kinship coefficient estimates in \code{kinMat} will also be used as divergence measures in place of \code{divMat}.
		
	It is important that the order of individuals in the matrices \code{kinMat} and \code{divMat} match the order of individuals in the \code{genoData}.  
	
	There are multiple ways to partition the sample into an ancestry representative 'unrelated subset' and a 'related subset'.  If \code{kinMat} is specified and \code{unrel.set = NULL}, then the PC-AiR algorithm is used to find an 'optimal' partition (see 'References' for a paper describing the algorithm).  If \code{kinMat = NULL} and \code{unrel.set} is specified, then the individuals with IDs in \code{unrel.set} are used as the 'unrelated subset'.  If both \code{kinMat} and \code{unrel.set} are specified, then all individuals with IDs in \code{unrel.set} are forced in the 'unrelated subset' and the PC-AiR algorithm is used to partition the rest of the sample; this is especially useful for including reference samples of known ancestry in the 'unrelated subset'.  If \code{kinMat = NULL} and \code{unrel.set = NULL}, then a standard principal components analysis that does not account for relatedness is performed.
}
\value{An object of class '\code{pcair}'.  A list including:
	\item{vectors}{A matrix of the top \code{v} principal components; each column is a principal component. Sample IDs are provided as rownames.}
	\item{values}{A vector of eigenvalues matching the top \code{v} principal components.  These values are determined from the standard PCA run on the 'unrelated subset'.}
	\item{sum.values}{The sum of all the eigenvalues from the standard PCA run on the 'unrelated subset' (regardless of how many were returned).}
	\item{rels}{A vector of IDs for individuals in the 'related subset'.}
	\item{unrels}{A vector of IDs for individuals in the 'unrelated subset'.}
	\item{kin.thresh}{The threshold value used for declaring each pair of individuals as related or unrelated.}
	\item{div.thresh}{The threshold value used for determining if each pair of individuals is ancestrally divergent.}
	\item{nsamp}{The total number of samples in the analysis.}
	\item{nsnps}{The total number of SNPs used in the analysis, after filtering on \code{MAF}.}
	\item{MAF}{The minor allele frequency (MAF) filter used on SNPs.}
	\item{call}{The function call passed to \code{pcair}.}
	\item{method}{A character string.  Either "PC-AiR" or "Standard PCA" identifying which method was used for computing principal components.}
}
\references{
	Conomos M.P., Miller M., & Thornton T. (2015). Robust Inference of Population Structure for Ancestry Prediction and Correction of Stratification in the Presence of Relatedness. Genetic Epidemiology, 39(4), 276-293.
    
    Gogarten, S.M., Bhangale, T., Conomos, M.P., Laurie, C.A., McHugh, C.P., Painter, I., ... & Laurie, C.C. (2012). GWASTools: an R/Bioconductor package for quality control and analysis of Genome-Wide Association Studies. Bioinformatics, 28(24), 3329-3331.
    
    Manichaikul, A., Mychaleckyj, J.C., Rich, S.S., Daly, K., Sale, M., & Chen, W.M. (2010). Robust relationship inference in genome-wide association studies. Bioinformatics, 26(22), 2867-2873.
}
\author{Matthew P. Conomos}
\note{The \code{GenotypeData} function in the \code{GWASTools} package should be used to create the input \code{genoData}.  Input to the \code{GenotypeData} function can easily be created from an R matrix or GDS file.  PLINK .bed, .bim, and .fam files can easily be converted to a GDS file with the function \code{snpgdsBED2GDS} in the \code{SNPRelate} package.}
\seealso{
	\code{\link{pcairPartition}} for a description of the function used by \code{pcair} that can be used to partition the sample into 'unrelated' and 'related' subsets without performing PCA.
	\code{\link{plot.pcair}} for plotting.
	\code{\link{king2mat}} for creating a matrix of pairwise kinship coefficient estimates from KING output text files that can be used for \code{kinMat} or \code{divMat}.
	\code{\link[GWASTools:GWASTools-package]{GWASTools}} for a description of the package containing the following functions: \code{\link[GWASTools]{GenotypeData}} for a description of creating a \code{GenotypeData} class object for storing sample and SNP genotype data, \code{\link[GWASTools]{MatrixGenotypeReader}} for a description of reading in genotype data stored as a matrix, and \code{\link[GWASTools]{GdsGenotypeReader}} for a description of reading in genotype data stored as a GDS file.  Also see \code{\link[SNPRelate]{snpgdsBED2GDS}} in the \code{\link[SNPRelate:SNPRelate-package]{SNPRelate}} package for a description of converting binary PLINK files to GDS.  The generic functions \code{\link{summary}} and \code{\link{print}}.
}
\examples{
# file path to GDS file
gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
# read in GDS data
HapMap_geno <- GdsGenotypeReader(filename = gdsfile)
# create a GenotypeData class object
HapMap_genoData <- GenotypeData(HapMap_geno)
# load saved matrix of KING-robust estimates
data("HapMap_ASW_MXL_KINGmat")
# run PC-AiR
mypcair <- pcair(genoData = HapMap_genoData, kinMat = HapMap_ASW_MXL_KINGmat, 
                divMat = HapMap_ASW_MXL_KINGmat)
close(HapMap_genoData)
}
\keyword{robust}
\keyword{multivariate}
\keyword{ancestry}

