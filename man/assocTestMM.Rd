\name{assocTestMM}
\alias{assocTestMM}
\title{SNP Genotype Association Testing with Mixed Models}
\description{\code{assocTestMM} performs SNP genotype association tests using the null model fit with \code{\link{fitNullMM}}.}
\usage{
assocTestMM(genoData, nullMMobj, test = "Wald", snp.include = NULL,
            chromosome = NULL, impute.geno = TRUE, snp.block.size = 5000,
            ivars = NULL, ivar.return.betaCov = FALSE, verbose = TRUE)
}
\arguments{
    \item{genoData}{An object of class \code{GenotypeData} from the package \code{GWASTools} containing the genotype data for SNPs and samples to be used for the analysis. This object can easily be created from a matrix of SNP genotype data, PLINK files, or GDS files.}
    \item{nullMMobj}{A null model object returned by \code{fitNullMM}.}
    \item{test}{A character string specifying the type of test to be performed. The possibilities are "Wald" (default) or "Score"; only "Score" can be used when the family of the null model fit with \code{fitNullMM} is not gaussian.}
    \item{snp.include}{A vector of SNP IDs to include in the analysis.  If NULL, see \code{chromosome} for further details.}
    \item{chromosome}{A vector of integers specifying which chromosomes to analyze.  This parameter is only considred when \code{snp.include} is NULL; if \code{chromosome} is also NULL, then all SNPs are included.}
    \item{impute.geno}{A logical indicator of whether sporadic missing genotype values should be mean imputed. The default is TRUE. See 'Details' for further information.}
    \item{snp.block.size}{The number of SNPs to read-in/analyze at once.  The default value is 5000. See 'Details' for further information regarding how this parameter works when \code{impute.geno} is FALSE.}
    \item{ivars}{A vector of character strings specifying the names of the variables for which a genotype interaction term should be included. If NULL (default) no genotype interactions are included. See 'Details' for further information.}
    \item{ivar.return.betaCov}{Logical indicator of whether the estimated covariance matrix of the effect size estimates (betas) for the genotype and genotype interaction terms should be returned; the default is FALSE.}
    \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is TRUE.}
}
\details{
    When \code{impute.geno} is TRUE, sporadic missing genotype values are mean imputed using the minor allele frequency (MAF) calculated on all other samples at that SNP. When \code{impute.geno} is FALSE, samples with missing values for all of the SNP genotypes in the current SNP block are removed from the analysis for the block; this may significantly slow down computation time because many pre-computed matrices need to be re-computed each time the sample set changes.  Also note: when \code{impute.geno} is FALSE, sporadic missingness for a sample inside of a SNP block will lead to an error.
    
    The input \code{ivars} can be used to perform GxE tests.  Multiple interaction variables may be specified, but all interaction variables specified must have been included as covariates in fitting the null model with \code{fitNullMM}.  When performing GxE analyses, \code{assocTestMM} will report two tests:  (1) the joint test of all genotype interaction terms in the model (this is the test for any genotype interaction effect), and (2) the joint test of the genotype term along with all of the genotype interaction terms (this is the test for any genetic effect). Individual genotype interaction terms can be tested by creating Wald test statistics from the reported effect size estimates and their standard errors (Note: when \code{ivars} contains a single continuous or binary covariate, this test is the same as the test for any genotype interaction effect mentioned above).  In order to test more complex hypotheses regarding subsets of multiple genotype interaction terms, \code{ivar.return.betaCov} can be used to retrieve the estimated covariance matrix of the effect size estimates.
}
\value{A data.frame where each row refers to a different SNP with the columns:
    \item{snpID}{The SNP ID}
    \item{chr}{The numeric chromosome value}
    \item{n}{The number of samples used to analyze the SNP}
    \item{MAF}{The estimated minor allele frequency}
    \item{minor.allele}{Either "A" or "B" indicating which allele is the minor allele}
    If \code{test} is "Score":
    \item{Score}{The value of the score function}
    \item{Var}{The variance of the score function}
    \item{Score.Stat}{The score chi-squared test statistic}
    \item{Score.pval}{The score p-value}
    If \code{test} is "Wald" and \code{ivars} is \code{NULL}:
    \item{Est}{The effect size estimate for each additional copy of the "A" allele}
    \item{SE}{The estimated standard error of the effect size estimate}    
    \item{Wald.Stat}{The Wald chi-squared test statistic}
    \item{Wald.pval}{The Wald p-value}
    If \code{test} is "Wald" and \code{ivars} is not \code{NULL}:
    \item{Est.G}{The effect size estimate for the genotype term}
    \item{Est.G:ivar}{The effect size estimate for the genotype*ivar interaction term. There will be as many of these terms as there are interaction variables, and "ivar" will be replaced with the variable name.}
    \item{SE.G}{The estimated standard error of the genotype term effect size estimate}
    \item{SE.G:ivar}{The estimated standard error of the genotype*ivar effect size estimate. There will be as many of these terms as there are interaction variables, and "ivar" will be replaced with the variable name.}
    \item{GxE.Stat}{The Wald chi-squared test statistic for the test of all genotype interaction terms.  When there is only one genotype interaction term, this is the test statistic for that term.}
    \item{GxE.pval}{The Wald p-value for the test of all genotype interaction terms; i.e. the test of any genotype interaction effect}
    \item{Joint.Stat}{The Wald chi-squared test statistic for the joint test of the genotype term and all of the genotype interaction terms}
    \item{Joint.pval}{The Wald p-value for the joint test of the genotype term and all of the genotype interaction terms; i.e. the test of any genotype effect}
    
    When \code{ivars} is not \code{NULL}, if \code{ivar.return.betaCov} is \code{TRUE}, then the output is a list with two elements.  The first, "results", is the data.frame described above.  The second, "betaCov", is a list with length equal to the number of rows of "results", where each element of the list is the covariance matrix of the effect size estimates (betas) for the genotype and genotype interaction terms.
}
%\references{}
\author{Matthew P. Conomos}
\note{The \code{GenotypeData} function in the \code{GWASTools} package should be used to create the input \code{genoData}.  Input to the \code{GenotypeData} function can easily be created from an R matrix or GDS file.  PLINK .bed, .bim, and .fam files can easily be converted to a GDS file with the function \code{snpgdsBED2GDS} in the \code{SNPRelate} package.}
\seealso{
     \code{\link{fitNullMM}} for fitting the null mixed model needed as input to \code{assocTestMM}.
     \code{\link[GWASTools]{qqPlot}} for a function to make QQ plots and \code{\link[GWASTools]{manhattanPlot}} for a function to make Manhattan plots of p-values.
     \code{\link[GWASTools:GWASTools-package]{GWASTools}} for a description of the package containing the following functions: \code{\link[GWASTools]{GenotypeData}} for a description of creating a \code{GenotypeData} class object for storing sample and SNP genotype data, \code{\link[GWASTools]{MatrixGenotypeReader}} for a description of reading in genotype data stored as a matrix, and \code{\link[GWASTools]{GdsGenotypeReader}} for a description of reading in genotype data stored as a GDS file.  Also see \code{\link[SNPRelate]{snpgdsBED2GDS}} in the \code{\link[SNPRelate:SNPRelate-package]{SNPRelate}} package for a description of converting binary PLINK files to GDS.
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
                
# run PC-Relate
mypcrel <- pcrelate(genoData = HapMap_genoData, pcMat = mypcair$vectors[,1],
        		training.set = mypcair$unrels)

# generate a phenotype
set.seed(4)
pheno <- 0.2*mypcair$vectors[,1] + rnorm(mypcair$nsamp, mean = 0, sd = 1)

# make ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID = mypcrel$sample.id, 
              pc1 = mypcair$vectors[,1], pheno = pheno))

# make covMatList
covMatList <- list("Kin" = pcrelateMakeGRM(mypcrel))

# fit the null mixed model
nullmod <- fitNullMM(scanData = scanAnnot, outcome = "pheno", covars = "pc1", covMatList = covMatList)

# run the association test
myassoc <- assocTestMM(genoData = HapMap_genoData, nullMMobj = nullmod)
close(HapMap_genoData)

# make a QQ plot
qqPlot(myassoc$Wald.pval)
}
\keyword{mixed model}
\keyword{association}
