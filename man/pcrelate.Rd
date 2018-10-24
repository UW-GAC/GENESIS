\name{pcrelate}
\alias{pcrelate}
\alias{pcrelate,GenotypeIterator-method}
\alias{pcrelate,SeqVarIterator-method}
\title{PC-Relate: Model-Free Estimation of Recent Genetic Relatedness}
\description{\code{pcrelate} is used to estimate kinship coefficients, IBD sharing probabilities, and inbreeding coefficients using genome-wide SNP data.  PC-Relate accounts for population structure (ancestry) among sample individuals through the use of ancestry representative principal components (PCs) to provide accurate relatedness estimates due only to recent family (pedigree) structure.}
\usage{
\S4method{pcrelate}{GenotypeIterator}(gdsobj,
	pcs,
	scale = c('overall', 'variant', 'none'),
	ibd.probs = TRUE,
	sample.include = NULL,
	training.set = NULL,
	sample.block.size = 5000,
	maf.thresh = 0.01,
	maf.bound.method = c('filter', 'truncate'),
	small.samp.correct = FALSE,
	num.cores = 1,
	verbose = TRUE)
\S4method{pcrelate}{SeqVarIterator}(gdsobj,
	pcs,
	scale = c('overall', 'variant', 'none'),
	ibd.probs = TRUE,
	sample.include = NULL,
	training.set = NULL,
	sample.block.size = 5000,
	maf.thresh = 0.01,
	maf.bound.method = c('filter', 'truncate'),
	small.samp.correct = FALSE,
	num.cores = 1,
	verbose = TRUE)
}
\arguments{
    \item{gdsobj}{An object of class \code{\link{SeqVarIterator}} from the package \pkg{\link{SeqVarTools}}, or an object of class \code{\link{GenotypeIterator}} from the package \pkg{\link{GWASTools}}, containing the genotype data for the variants and samples to be used for the analysis.}
  \item{pcs}{A matrix of principal components (PCs) to be used for ancestry adjustment.  Each column represents a PC, and each row represents an individual. IDs for each individual must be set as the row names of the matrix.}
  \item{scale}{A character string taking the values 'overall', 'variant', or 'none' indicating how genotype values should be standardized.  This should be set to 'overall' (the default) in order to do a PC-Relate analysis; see 'Details' for more information.}
  \item{ibd.probs}{Logical indicator of whether pairwise IBD sharing probabilities (k0, k1, k2) should be estimated; the default is TRUE.}
  \item{sample.include}{A vector of IDs for samples to include in the analysis.  If NULL, all samples in \code{gdsobj} are included.}
  \item{training.set}{An optional vector of IDs identifying which samples to use for estimation of the ancestry effect when estimating individual-specific allele frequencies.  If NULL, all samples in sample.include are used. See 'Details' for more information.}
  \item{sample.block.size}{The number of individuals to read-in/analyze at once; the default value is 5000. See 'Details' for more information.}
  \item{maf.thresh}{Minor allele frequency threshold; if an individual's estimated individual-specific minor allele frequency at a SNP is less than this value, that indivdiual will either have that SNP excluded from the analysis or have their estimated indivdiual-specific minor allele frequency truncated to this value, depending on \code{maf.bound.method}. The default value is 0.01.}
  \item{maf.bound.method}{How individual-specific minor allele frequency estimates less that \code{maf.thresh} are handled. When set to `filter` (default), SNPs for which an individual's estimated individual-specific minor allele frequency are below \code{maf.thresh} are excluded from the analysis for that individual. When set to `truncate`, estimated individual-specific minor allele frequncies below \code{maf.thresh} have their value set to \code{maf.thresh}.}
 \item{small.samp.correct}{Logical indicator of whether to implement a small sample correction.}
  \item{num.cores}{The number of cores to use. If \code{> 1}, multiple blocks of variants and/or samples can be analyzed in parallel.}
  \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is TRUE.}
 }
\details{
	The basic premise of PC-Relate is to estimate kinship coefficients, IBD sharing probabilities, and inbreeding coefficients that reflect recent family (pedigree) relatedness by conditioning out genetic similarity due to distant population structure (ancestry) with ancestry representative principal components (PCs).

	It is important that the PCs used in \code{pcs} to adjust for ancestry are representative of ancestry and NOT family structure, so we recommend using PCs calculated with PC-AiR (see: \code{\link{pcair}}).

    In order to perform relatedness estimation, allele frequency estimates are required for centering and scaling genotype values. Individual-specific allele frequencies calculated for each individual at each SNP using the PCs specified in \code{pcs} are used. There are muliple choices for how genotype values are scaled. When \code{scale} is 'variant', centered genotype values at each SNP are divided by their expected variance under Hardy-Weinberg equilibrium. When \code{scale} is 'overall', centered genotype values at all SNPs are divided by the average across all SNPs of their expected variances under Hardy-Weinberg equilibrium; this scaling leads to more stable behavior when using low frequency variants. When \code{scale} is 'none', genotype values are only centered and not scaled; this won't provide accurate kinship coefficient estimates but may be useful for other purposes. Set \code{scale} to 'overall' to perform a standard PC-Relate analysis; this is the default.  If \code{scale} is set to 'variant', the estimators are very similar to REAP.

	The optional input \code{training.set} allows the user to specify which samples are used to estimate the ancestry effect when estimating individual-specific allele frequencies. Ideally, \code{training.set} is a set of mutually unrelated individuals.  If prior information regarding pedigree structure is available, this can be used to select \code{training.set}, or if \code{\link{pcair}} was used to obtain the PCs, then the individuals in the PC-AiR 'unrelated subset' can be used.  If no prior information is available, all individuals should be used.

	The \code{sample.block.size} can be specified to alleviate memory issues when working with very large data sets.  If \code{sample.block.size} is smaller than the number of individuals included in the analysis, then individuals will be analyzed in separate blocks.  This reduces the memory required for the analysis, but genotype data must be read in multiple times for each block (to analyze all pairs), which increases the number of computations required.
}
\value{An object of class '\code{pcrelate}'. A list including:
	\item{kinBtwn}{A data.frame of estimated pairwise kinship coefficients and IBD sharing probabilities (if \code{ibd.probs} is \code{TRUE}).}
	\item{kinSelf}{A data.frame of estimated inbreeding coefficients.}
	%\item{ibd.probs}{A matrix of estimated pairwise IBD sharing probabilities; the lower triangle gives k0 (the probability of sharing 0 alleles IBD), the upper triangle gives k2 (the probability of sharing 2 alleles IBD), and the diagonal is missing. The order of samples matches \code{sample.id}. This matrix is returned only if \code{ibd.probs = TRUE} in the input.}
	%\item{nsnp}{A matrix specifying the the number of SNPs used to estimate the relatedness measures for each pair of individuals. The order of samples matches \code{sample.id}.}
	%\item{kincorrect}{A vector specifying the correction factors used for the small sample correction, or NULL.}
	%\item{k2correct}{A vector specifying the correction factors used for the small sample correction, or NULL.}
}
\references{
	Conomos M.P., Reiner A.P., Weir B.S., & Thornton T.A. (2016). Model-free Estimation of Recent Genetic Relatedness. American Journal of Human Genetics, 98(1), 127-148.
}
\author{Matthew P. Conomos}
\seealso{
	\code{\link{pcrelateToMatrix}}
}
\examples{
library(GWASTools)

# file path to GDS file
gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
# read in GDS data
HapMap_geno <- GdsGenotypeReader(filename = gdsfile)
# create a GenotypeData class object
HapMap_genoData <- GenotypeData(HapMap_geno)
# load saved matrix of KING-robust estimates
data("HapMap_ASW_MXL_KINGmat")
# run PC-AiR
mypcair <- pcair(HapMap_genoData, kinobj = HapMap_ASW_MXL_KINGmat, 
                divobj = HapMap_ASW_MXL_KINGmat)

# create a GenotypeBlockIterator object
HapMap_genoData <- GenotypeBlockIterator(HapMap_genoData)
# run PC-Relate
mypcrel <- pcrelate(HapMap_genoData, pcs = mypcair$vectors[,1,drop=FALSE],
				training.set = mypcair$unrels)
head(mypcrel$kinBwtn)
head(mypcrel$kinSelf)

grm <- pcrelateToMatrix(mypcrel)
dim(grm)

close(HapMap_genoData)
}
\keyword{robust}
\keyword{relatedness}
