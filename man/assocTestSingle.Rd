\name{assocTestSingle}
\alias{assocTestSingle}
\alias{assocTestSingle-methods}
\alias{assocTestSingle,SeqVarIterator-method}

\title{Genotype Association Testing with Mixed Models}

\description{\code{assocTestSingle} performs genotype association tests
  using the null model fit with \code{\link{fitNullModel}}.}

\usage{
\S4method{assocTestSingle}{SeqVarIterator}(gdsobj, null.model, test=c("Score", "Wald"),
                GxE=NULL, verbose=TRUE)
}

\arguments{
    \item{gdsobj}{An object of class \code{\link{SeqVarIterator}} from the package \pkg{\link{SeqVarTools}} containing the genotype data for the variants and samples to be used for the analysis.}
    \item{null.model}{A null model object returned by \code{\link{fitNullModel}}.}
    \item{test}{A character string specifying the type of test to be performed. The possibilities are \code{"Score"} (default) or \code{"Wald"}; only \code{"Score"} can be used when the family of the null model fit with \code{\link{fitNullModel}} is not gaussian.}
    %\item{impute.geno}{A logical indicator of whether sporadic missing genotype values should be mean imputed. The default is TRUE. See 'Details' for further information.}
    \item{GxE}{A vector of character strings specifying the names of the variables for which a genotype interaction term should be included. If \code{NULL} (default) no genotype interactions are included. See 'Details' for further information.}
    %\item{ivar.return.betaCov}{Logical indicator of whether the estimated covariance matrix of the effect size estimates (betas) for the genotype and genotype interaction terms should be returned; the default is FALSE.}
    \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is \code{TRUE}.}
}

\details{
    %When \code{impute.geno} is TRUE, sporadic missing genotype values are mean imputed using the minor allele frequency (MAF) calculated on all other samples at that SNP. When \code{impute.geno} is FALSE, samples with missing values for all of the SNP genotypes in the current SNP block are removed from the analysis for the block; this may significantly slow down computation time because many pre-computed matrices need to be re-computed each time the sample set changes.  Also note: when \code{impute.geno} is FALSE, sporadic missingness for a sample inside of a SNP block will lead to an error.

    Sporadic missing genotype values are mean imputed using the minor allele frequency (MAF) calculated on all other samples at that variant.
    
    The input \code{GxE} can be used to perform GxE tests.  Multiple interaction variables may be specified, but all interaction variables specified must have been included as covariates in fitting the null model with \code{fitNullModel}.  When performing GxE analyses, \code{assocTestSingle} will report two tests:  (1) the joint test of all genotype interaction terms in the model (this is the test for any genotype interaction effect), and (2) the joint test of the genotype term along with all of the genotype interaction terms (this is the test for any genetic effect). Individual genotype interaction terms can be tested by creating Wald test statistics from the reported effect size estimates and their standard errors (Note: when \code{GxE} contains a single continuous or binary covariate, this test is the same as the test for any genotype interaction effect mentioned above).  %In order to test more complex hypotheses regarding subsets of multiple genotype interaction terms, \code{ivar.return.betaCov} can be used to retrieve the estimated covariance matrix of the effect size estimates.
}

\value{A data.frame where each row refers to a different variant with the columns:
    \item{variant.id}{The variant ID}
    \item{chr}{The chromosome value}
    \item{pos}{The base pair position}
    \item{allele.index}{The index of the alternate allele. For biallelic variants, this will always be 1.}
    \item{n.obs}{The number of samples with non-missing genotypes}
    \item{freq}{The estimated alternate allele frequency}
    If \code{test} is \code{"Score"}:
    \item{Score}{The value of the score function}
    \item{Score.SE}{The estimated standard error of the Score}
    \item{Score.Stat}{The score Z test statistic}
    \item{Score.pval}{The score p-value}
    If \code{test} is \code{"Wald"} and \code{GxE} is \code{NULL}:
    \item{Est}{The effect size estimate for each additional copy of the alternate allele}
    \item{Est.SE}{The estimated standard error of the effect size estimate}    
    \item{Wald.Stat}{The Wald Z test statistic}
    \item{Wald.pval}{The Wald p-value}
    If \code{test} is \code{"Wald"} and \code{GxE} is not \code{NULL}:
    \item{Est.G}{The effect size estimate for the genotype term}
    \item{Est.G:env}{The effect size estimate for the genotype*env interaction term. There will be as many of these terms as there are interaction variables, and "env" will be replaced with the variable name.}
    \item{SE.G}{The estimated standard error of the genotype term effect size estimate}
    \item{SE.G:env}{The estimated standard error of the genotype*env effect size estimate. There will be as many of these terms as there are interaction variables, and "env" will be replaced with the variable name.}
    \item{GxE.Stat}{The Wald Z test statistic for the test of all genotype interaction terms.  When there is only one genotype interaction term, this is the test statistic for that term.}
    \item{GxE.pval}{The Wald p-value for the test of all genotype interaction terms; i.e. the test of any genotype interaction effect}
    \item{Joint.Stat}{The Wald Z test statistic for the joint test of the genotype term and all of the genotype interaction terms}
    \item{Joint.pval}{The Wald p-value for the joint test of the genotype term and all of the genotype interaction terms; i.e. the test of any genotype effect}
    
    %When \code{GxE} is not \code{NULL}, if \code{ivar.return.betaCov} is \code{TRUE}, then the output is a list with two elements.  The first, "results", is the data.frame described above.  The second, "betaCov", is a list with length equal to the number of rows of "results", where each element of the list is the covariance matrix of the effect size estimates (betas) for the genotype and genotype interaction terms.

    The effect size estimate is for each copy of the alternate allele.
    For multiallelic variants, each alternate allele is tested separately.
}

\author{Matthew P. Conomos, Stephanie M. Gogarten, Tamar Sofer, Ken Rice, Chaoyu Yu}

\seealso{
     \code{\link{fitNullModel}} for fitting the null mixed model needed as input to \code{assocTestSingle}.
     \code{\link{SeqVarIterator}} for creating the input object with genotypes.
     %\code{\link[GWASTools]{qqPlot}} for a function to make QQ plots and \code{\link[GWASTools]{manhattanPlot}} for a function to make Manhattan plots of p-values.
     %\code{\link[GWASTools:GWASTools-package]{GWASTools}} for a description of the package containing the following functions: \code{\link{GenotypeData}} for a description of creating a \code{GenotypeData} class object for storing sample and SNP genotype data, \code{\link{MatrixGenotypeReader}} for a description of reading in genotype data stored as a matrix, and \code{\link{GdsGenotypeReader}} for a description of reading in genotype data stored as a GDS file.  Also see \code{\link[SNPRelate]{snpgdsBED2GDS}} in the \code{\link[SNPRelate:SNPRelate-package]{SNPRelate}} package for a description of converting binary PLINK files to GDS.
}

\examples{
library(SeqVarTools)
library(Biobase)

# open a sequencing GDS file
gdsfile <- seqExampleFileName("gds")
gds <- seqOpen(gdsfile)

# simulate some phenotype data
data(pedigree)
pedigree <- pedigree[match(seqGetData(gds, "sample.id"), pedigree$sample.id),]
pedigree$outcome <- rnorm(nrow(pedigree))

# construct a SeqVarIterator object
seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(pedigree))
iterator <- SeqVarBlockIterator(seqData)

# fit the null model
nullmod <- fitNullModel(iterator, outcome="outcome", covars="sex")

# run the association test
assoc <- assocTestSingle(iterator, nullmod)

seqClose(iterator)
}
\keyword{mixed model}
\keyword{association}
