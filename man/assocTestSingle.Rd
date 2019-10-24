\name{assocTestSingle}
\alias{assocTestSingle}
\alias{assocTestSingle-methods}
\alias{assocTestSingle,SeqVarIterator-method}
\alias{assocTestSingle,GenotypeIterator-method}

\title{Genotype Association Testing with Mixed Models}

\description{\code{assocTestSingle} performs genotype association tests
  using the null model fit with \code{\link{fitNullModel}}.}

\usage{
\S4method{assocTestSingle}{SeqVarIterator}(gdsobj, null.model, test=c("Score", "Score.SPA"),
                recalc.pval.thresh=0.05, GxE=NULL, sparse=TRUE, imputed=FALSE, male.diploid=TRUE,
                genome.build=c("hg19", "hg38"), verbose=TRUE)
\S4method{assocTestSingle}{GenotypeIterator}(gdsobj, null.model, test=c("Score", "Score.SPA"),
                recalc.pval.thresh=0.05, GxE=NULL, male.diploid=TRUE, verbose=TRUE)
}

\arguments{
    \item{gdsobj}{An object of class \code{\link{SeqVarIterator}} from the package \pkg{\link{SeqVarTools}}, or an object of class \code{\link{GenotypeIterator}} from the package \pkg{\link{GWASTools}}, containing the genotype data for the variants and samples to be used for the analysis.}
    \item{null.model}{A null model object returned by \code{\link{fitNullModel}}.}
    \item{test}{A character string specifying the type of test to be performed. The possibilities are \code{"Score"} (default) and \code{"Score.SPA"}; \code{"Score.SPA"} can only be used when the family of the null model fit with \code{\link{fitNullModel}} is binomial.}
    %\item{impute.geno}{A logical indicator of whether sporadic missing genotype values should be mean imputed. The default is TRUE. See 'Details' for further information.}
    \item{recalc.pval.thresh}{If test is not "Score",  recalculate p-values using the specified `test` for variants with a Score p-value below this threshold; return the score p-value for all other variants.}
    \item{GxE}{A vector of character strings specifying the names of the variables for which a genotype interaction term should be included. If \code{NULL} (default) no genotype interactions are included. See 'Details' for further information.}
    %\item{ivar.return.betaCov}{Logical indicator of whether the estimated covariance matrix of the effect size estimates (betas) for the genotype and genotype interaction terms should be returned; the default is FALSE.}
    \item{sparse}{Logical indicator of whether to read genotypes as sparse Matrix objects; the default is \code{TRUE}. Set this to \code{FALSE} if the alternate allele dosage of the genotypes in the test are not expected to be mostly 0.}
    \item{imputed}{Logical indicator of whether to read dosages from the "DS" field containing imputed dosages instead of counting the number of alternate alleles.}
    \item{male.diploid}{Logical for whether males on sex chromosomes are coded as diploid.}
    \item{genome.build}{A character sting indicating genome build; used to identify pseudoautosomal regions on the X and Y chromosomes.}
    \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is \code{TRUE}.}
}

\details{
  All samples included in \code{null model} will be included in the association test, even if a different set of samples is present in the current filter for \code{gdsobj}.
  
    The effect size estimate is for each copy of the alternate allele. For multiallelic variants, each alternate allele is tested separately.
    %When \code{impute.geno} is TRUE, sporadic missing genotype values are mean imputed using the minor allele frequency (MAF) calculated on all other samples at that SNP. When \code{impute.geno} is FALSE, samples with missing values for all of the SNP genotypes in the current SNP block are removed from the analysis for the block; this may significantly slow down computation time because many pre-computed matrices need to be re-computed each time the sample set changes.  Also note: when \code{impute.geno} is FALSE, sporadic missingness for a sample inside of a SNP block will lead to an error.

    Sporadic missing genotype values are mean imputed using the allele frequency calculated on all other samples at that variant.

    Monomorphic variants (including variants where every sample is a heterozygote) are omitted from the results.
 
    The input \code{GxE} can be used to perform GxE tests.  Multiple interaction variables may be specified, but all interaction variables specified must have been included as covariates in fitting the null model with \code{fitNullModel}.  When performing GxE analyses, \code{assocTestSingle} will report two tests:  (1) the joint test of all genotype interaction terms in the model (this is the test for any genotype interaction effect), and (2) the joint test of the genotype term along with all of the genotype interaction terms (this is the test for any genetic effect). Individual genotype interaction terms can be tested by creating test statistics from the reported effect size estimates and their standard errors (Note: when \code{GxE} contains a single continuous or binary covariate, this test is the same as the test for any genotype interaction effect mentioned above).  %In order to test more complex hypotheses regarding subsets of multiple genotype interaction terms, \code{ivar.return.betaCov} can be used to retrieve the estimated covariance matrix of the effect size estimates.

    The saddle point approximation (SPA), run by using \code{test = "Score.SPA"}, implements the method described by Dey et al. (2017), which was extended to mixed models by Zhou et al. (2018) in their SAIGE software. SPA provides better calibration of p-values when fitting mixed models with a binomial family for a sample with an imbalanced case to control ratio.

    For the \code{\link{GenotypeIterator}} method, objects created with \code{\link{GdsGenotypeReader}} or \code{\link{MatrixGenotypeReader}} are supported. \code{\link{NcdfGenotypeReader}} objects are not supported.
}

\value{A data.frame where each row refers to a different variant with the columns:
    \item{variant.id}{The variant ID}
    \item{chr}{The chromosome value}
    \item{pos}{The base pair position}
    \item{allele.index}{The index of the alternate allele. For biallelic variants, this will always be 1.}
    \item{n.obs}{The number of samples with non-missing genotypes}
    \item{freq}{The estimated alternate allele frequency}
    \item{MAC}{The minor allele count. For multiallelic variants, "minor" is determined by comparing the count of the alternate allele specified by \code{allele.index} with the sum of all other alleles.}
    If \code{test} is \code{"Score"}:
    \item{Score}{The value of the score function}
    \item{Score.SE}{The estimated standard error of the Score}
    \item{Score.Stat}{The score Z test statistic}
    \item{Score.pval}{The score p-value}
    \item{Est}{An approximation of the effect size estimate for each additional copy of the alternate allele}
    \item{Est.SE}{An approximation of the standard error of the effect size estimate}
    \item{PVE}{An approximation of the proportion of phenotype variance explained}
    % If \code{test} is \code{"Wald"} and \code{GxE} is \code{NULL}:
    % \item{Est}{The effect size estimate for each additional copy of the alternate allele}
    % \item{Est.SE}{The estimated standard error of the effect size estimate}    
    % \item{Wald.Stat}{The Wald Z test statistic}
    % \item{Wald.pval}{The Wald p-value}
    If \code{test} is \code{"Score.SPA"}:
    \item{SPA.pval}{The score p-value after applying the saddle point approximation (SPA)}
    \item{SPA.converged}{logical indiactor of whether the SPA converged; \code{NA} indicates that the SPA was not applied and the original Score.pval was returned}
    If \code{GxE} is not \code{NULL}:
    \item{Est.G}{The effect size estimate for the genotype term}
    \item{Est.G:env}{The effect size estimate for the genotype*env interaction term. There will be as many of these terms as there are interaction variables, and "env" will be replaced with the variable name.}
    \item{SE.G}{The estimated standard error of the genotype term effect size estimate}
    \item{SE.G:env}{The estimated standard error of the genotype*env effect size estimate. There will be as many of these terms as there are interaction variables, and "env" will be replaced with the variable name.}
    \item{GxE.Stat}{The Wald Z test statistic for the test of all genotype interaction terms.  When there is only one genotype interaction term, this is the test statistic for that term.}
    \item{GxE.pval}{The Wald p-value for the test of all genotype interaction terms; i.e. the test of any genotype interaction effect}
    \item{Joint.Stat}{The Wald Z test statistic for the joint test of the genotype term and all of the genotype interaction terms}
    \item{Joint.pval}{The Wald p-value for the joint test of the genotype term and all of the genotype interaction terms; i.e. the test of any genotype effect}
    
    %When \code{GxE} is not \code{NULL}, if \code{ivar.return.betaCov} is \code{TRUE}, then the output is a list with two elements.  The first, "results", is the data.frame described above.  The second, "betaCov", is a list with length equal to the number of rows of "results", where each element of the list is the covariance matrix of the effect size estimates (betas) for the genotype and genotype interaction terms.
}

\references{
    Dey, R., Schmidt, E. M., Abecasis, G. R., & Lee, S. (2017). A fast and accurate algorithm to test for binary phenotypes and its application to PheWAS. The American Journal of Human Genetics, 101(1), 37-49.
    Zhou, W., Nielsen, J. B., Fritsche, L. G., Dey, R., Gabrielsen, M. E., Wolford, B. N., ... & Bastarache, L. A. (2018). Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies. Nature genetics, 50(9), 1335.
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
set.seed(4)
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


library(GWASTools)

# open a SNP-based GDS file
gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
gds <- GdsGenotypeReader(filename = gdsfile)

# simulate some phenotype data
set.seed(4)
pheno <- data.frame(scanID=getScanID(gds),
                    outcome=rnorm(nscan(gds)))

# construct a GenotypeIterator object
genoData <- GenotypeData(gds, scanAnnot=ScanAnnotationDataFrame(pheno))
iterator <- GenotypeBlockIterator(genoData)

# fit the null model
nullmod <- fitNullModel(iterator, outcome="outcome")

# run the association test
assoc <- assocTestSingle(iterator, nullmod)

close(iterator)
}
\keyword{mixed model}
\keyword{association}
