\name{pcrelate}
\alias{pcrelate}
\title{PC-Relate: Model-Free Estimation of Recent Genetic Relatedness}
\description{\code{pcrelate} is used to estimate kinship coefficients, IBD sharing probabilities, and inbreeding coefficients using genome-wide SNP data.  PC-Relate accounts for population structure (ancestry) among sample individuals through the use of ancestry representative principal components (PCs) to provide accurate relatedness estimates due only to recent family (pedigree) structure.}
\usage{
pcrelate(genoData, pcMat = NULL, freq.type = "individual", scale = "overall",
            ibd.probs = TRUE, scan.include = NULL, training.set = NULL, scan.block.size = 5000, 
			snp.include = NULL, chromosome = NULL, snp.block.size = 10000, 
			MAF = 0.01, write.to.gds = FALSE, gds.prefix = NULL, 
			correct = TRUE, verbose = TRUE)
}
\arguments{
  \item{genoData}{An object of class \code{GenotypeData} from the package \code{GWASTools} containing the genotype data for SNPs and samples to be used for the analysis.  This object can easily be created from a matrix of SNP genotype data, PLINK files, or GDS files.  Alternatively, this could be an object of class \code{SeqVarData} from the package \code{SeqVarTools} containing the genotype data for the sequencing variants and samples to be used for the analysis.} 
  \item{pcMat}{An optional matrix of principal components (PCs) to be used for ancestry adjustment.  Each column represents a PC, and each row represents an individual.  IDs for each individual must be set as the row names of the matrix.}
  \item{freq.type}{A character string taking the values 'individual' or 'population' indicating whether genotype values should be adjusted by individual-specific allele frequencies or population average allele frequencies. This should be set to 'individual' (the default) in order to do a PC-Relate analysis; see 'Details' for more information.}
  \item{scale}{A character string taking the values 'overall', 'variant', or 'none' indicating how genotype values should be standardized.  This should be set to 'overall' (the default) in order to do a PC-Relate analysis; see 'Details' for more information.}
  \item{ibd.probs}{Logical indicator of whether pairwise IBD sharing probabilities (k0, k1, k2) should be estimated; the default is TRUE.}
  \item{scan.include}{A vector of IDs for samples to include in the analysis.  If NULL, all samples in \code{genoData} are included.}
  \item{training.set}{An optional vector of IDs identifying which samples to use for estimation of the ancestry effect when estimating individual-specific allele frequencies.  If NULL, all samples in scan.include are used. See 'Details' for more information.}
  \item{scan.block.size}{The number of individuals to read-in/analyze at once; the default value is 5000. See 'Details' for more information.}
  \item{snp.include}{A vector of SNP IDs to include in the analysis.  If NULL, see \code{chromosome} for further details.}
  \item{chromosome}{A vector of integers specifying which chromosomes to analyze.  This parameter is only considred when \code{snp.include} is NULL; if \code{chromosome} is also NULL, then all SNPs are included.}
  \item{snp.block.size}{The number of SNPs to read-in/analyze at once.  The default value is 10000. }
  \item{MAF}{Minor allele frequency filter. When \code{freq.type} is 'individual', if an individual's estimated individual-specific minor allele frequency at a SNP is less than this value, that SNP will be excluded from the analysis for that individual. When \code{freq.type} is 'population', any SNP with a population minor allele frequency less than this value will be excluded from the analysis. The default value is 0.01.}
  \item{write.to.gds}{Logical indicator of whether the output should be written to GDS files.  If FALSE (the default), then the output is returned to the R console as expected.  See 'Details' for more information.}
  \item{gds.prefix}{File path specifying where to save the output when \code{write.to.gds = TRUE}.  If NULL, the prefix 'tmp' is used. See 'Details' for more information.}
  \item{correct}{Logical indicator of whether to implement a small sample correction.}
  \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is TRUE.}
 }
\details{
	The basic premise of PC-Relate is to estimate kinship coefficients, IBD sharing probabilities, and inbreeding coefficients that reflect recent family (pedigree) relatedness by conditioning out genetic similarity due to distant population structure (ancestry) with ancestry representative principal components (PCs).

	It is important that the PCs used in \code{pcMat} to adjust for ancestry are representative of ancestry and NOT family structure, so we recommend using PCs calculated with PC-AiR.

	It is important that the order of individuals in the matrix \code{pcMat} matches the order of individuals in \code{genoData}.
    
    In order to perform relatedness estimation, allele frequency estimates are required for centering and scaling genotype values. When \code{freq.type} is 'individual', individual-specific allele frequencies calculated for each individual at each SNP using the PCs specified in \code{pcMat} are used. When \code{freq.type} is 'population', population average allele frequencies calculated at each SNP are used for all individuals. (Note that when \code{freq.type} is set to 'population' there is no ancestry adjustment and the relatedness estimates will be confounded with population structure (ancestry)). There are muliple choices for how genotype values are scaled. When \code{scale} is 'variant', centered genotype values at each SNP are divided by their expected variance under Hardy-Weinberg equilibrium. When \code{scale} is 'overall', centered genotype values at all SNPs are divided by the average across all SNPs of their expected variances under Hardy-Weinberg equilibrium; this scaling leads to more stable behavior when using low frequency variants. When \code{scale} is 'none', genotype values are only centered and not scaled; this won't provide accurate kinship coefficient estimates but may be useful for other purposes. At a particular SNP, the variance used for scaling is either calculated separately for each individual using their individual-specific allele frequncies (when \code{freq.type} is 'individual') or once for all individuals using the population average allele frequency (when \code{freq.type} is 'population'). Set \code{freq.type} to 'individual' and \code{scale} to 'overall' to perform a standard PC-Relate analysis; these are the defaults.  If \code{freq.type} is set to 'individual' and \code{scale} is set to 'variant', the estimators are very similar to REAP. If \code{freq.type} is set to 'population' and \code{scale} is set to 'variant', the estimators are very similar to EIGENSOFT.

	The optional input \code{training.set} allows the user to specify which samples are used to estimate the ancestry effect when estimating individual-specific allele frequencies (if \code{freq.type} is 'individual') or to estimate the population allele frequency (if \code{freq.type} is 'population'. Ideally, \code{training.set} is a set of mutually unrelated individuals.  If prior information regarding pedigree structure is available, this can be used to select \code{training.set}, or if \code{\link{pcair}} was used to obtain the PCs, then the individuals in the PC-AiR 'unrelated subset' can be used.  If no prior information is available, all individuals should be used.

	The \code{scan.block.size} can be specified to alleviate memory issues when working with very large data sets.  If \code{scan.block.size} is smaller than the number of individuals included in the analysis, then individuals will be analyzed in separate blocks.  This reduces the memory required for the analysis, but genotype data must be read in multiple times for each block (to analyze all pairs), which increases the number of computations required.  NOTE: if individuals are broken up into more than 1 block, \code{write.to.gds} must be TRUE (see below).

	If \code{write.to.gds = TRUE}, then the output is written to two GDS files rather than returned to the R console. Use of this option requires the \code{\link[gdsfmt:gdsfmt-package]{gdsfmt}} package. The first GDS file, named ``<gds.prefix>_freq.gds'', contains the individual-specific allele frequency estimates for each individual at each SNP (when \code{freq.type} is 'individual') or the population allele frequency estimates at each SNP (when \code{freq.type} is 'population'.  The second GDS file, named ``<gds.prefix>_pcrelate.gds'', contains the PC-Relate output as described in Value below.
}
\value{An object of class '\code{pcrelate}'. A list including:
	\item{sample.id}{A vector of IDs for samples included in the analysis.}
	\item{kinship}{A matrix of estimated pairwise kinship coefficients. The order of samples matches \code{sample.id}.}
	\item{ibd.probs}{A matrix of estimated pairwise IBD sharing probabilities; the lower triangle gives k0 (the probability of sharing 0 alleles IBD), the upper triangle gives k2 (the probability of sharing 2 alleles IBD), and the diagonal is missing. The order of samples matches \code{sample.id}. This matrix is returned only if \code{ibd.probs = TRUE} in the input.}
	\item{nsnp}{A matrix specifying the the number of SNPs used to estimate the relatedness measures for each pair of individuals. The order of samples matches \code{sample.id}.}
	\item{kincorrect}{A vector specifying the correction factors used for the small sample correction, or NULL.}
	\item{k2correct}{A vector specifying the correction factors used for the small sample correction, or NULL.}
	\item{call}{The function call passed to \code{pcrelate}.}
	\item{method}{A character string.  Either 'PC-Relate' or 'Unadjusted' identifying which method was used for computing relatedness estimates.  'Unadjusted' is used when \code{pcMat = NULL} and corresponds to an assumption of population homogeneity.}
}
\references{
	Conomos M.P., Reiner A.P., Weir B.S., & Thornton T.A. (2016). Model-free Estimation of Recent Genetic Relatedness. American Journal of Human Genetics, 98(1), 127-148.
    
    Gogarten, S.M., Bhangale, T., Conomos, M.P., Laurie, C.A., McHugh, C.P., Painter, I., ... & Laurie, C.C. (2012). GWASTools: an R/Bioconductor package for quality control and analysis of Genome-Wide Association Studies. Bioinformatics, 28(24), 3329-3331.
}
\author{Matthew P. Conomos}
\note{The \code{GenotypeData} function in the \code{GWASTools} package should be used to create the input \code{genoData}.  Input to the \code{GenotypeData} function can easily be created from an R matrix or GDS file.  PLINK .bed, .bim, and .fam files can easily be converted to a GDS file with the function \code{snpgdsBED2GDS} in the \code{SNPRelate} package.  Alternatively, the \code{SeqVarData} function in the \code{SeqVarTools} package can be used to create the input \code{genodata} when working with sequencing data.}
\seealso{
	\code{\link{pcrelateReadKinship}}, \code{\link{pcrelateReadInbreed}}, and \code{\link{pcrelateMakeGRM}} for functions that can be used to read in the results output by \code{pcrelate}.
	\code{\link[GWASTools:GWASTools-package]{GWASTools}} for a description of the package containing the following functions: \code{\link[GWASTools]{GenotypeData}} for a description of creating a \code{GenotypeData} class object for storing sample and SNP genotype data, \code{\link[GWASTools]{MatrixGenotypeReader}} for a description of reading in genotype data stored as a matrix, and \code{\link[GWASTools]{GdsGenotypeReader}} for a description of reading in genotype data stored as a GDS file.  Also see \code{\link[SNPRelate]{snpgdsBED2GDS}} in the \code{\link[SNPRelate:SNPRelate-package]{SNPRelate}} package for a description of converting binary PLINK files to GDS.    
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
mypcair <- pcair(genoData = HapMap_genoData, kinMat = HapMap_ASW_MXL_KINGmat, 
                divMat = HapMap_ASW_MXL_KINGmat)
# run PC-Relate
mypcrel <- pcrelate(genoData = HapMap_genoData, pcMat = mypcair$vectors[,1],
				training.set = mypcair$unrels)
close(HapMap_genoData)
}
\keyword{robust}
\keyword{relatedness}
