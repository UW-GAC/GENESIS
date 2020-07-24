\name{jointScoreTest}
\alias{jointScoreTest}
\title{Perform a joint score test}
\description{\code{jointScoreTest} is used to perform a joint score test of a set of variants using a null model and a matrix of genotype dosages.}
\usage{
jointScoreTest(nullmod, G)
}
\arguments{
  \item{nullmod}{A null model object returned by \code{\link{fitNullModel}}.}
  \item{G}{A matrix of genotype dosages, where samples are the rows and variants are the columns.}
}
\details{
  \code{jointScoreTest} takes the object returned by \code{\link{fitNullModel}} and performs a joint score test for all variants in the \code{G} matrix using this null model.
  All effect size and PVE estimates in \code{fixef} are adjusted for (i.e. conditional on) all other variants included in \code{G}.

  The \code{G} matrix must be formatted such that the rows are samples and the columns are variants, and the entries are the dosage for that sample and variant.
  The matrix must have \code{sample.id} as rownames; this is used to match the genotypes to the null model.
  Therefore, the same sample identifiers must be used in both \code{nullmod} and \code{G}.
  \code{G} can contain additional samples and ordering is unimportant as long as all samples from the null model are present; it will be subset and reordered to match the set of samples in the null model.

  Colnames for \code{G} are not required.
  If present, the column names of \code{G} will be used as the rownames of the \code{fixef} element and the column and rownames of the \code{betaCov} element of the output.
  \code{fixef} and \code{betaCov} will be in the same order as the columns in \code{G}.

  Missing data in \code{G} are not allowed.
}
\value{\code{jointScoreTest} returns a list with the following elements:
  \item{Stat.joint}{Squared joint score test statistic for all variants in \code{G}.}
  \item{pval.joint}{p-value associated with the joint score test statistic drawn from a \deqn{\chi^2} distribution with \eqn{n_\mathrm{variants}} degrees of freedom.}
  \item{pve.joint}{Estimate of the proportion of variance explained by the variants in \code{G}.}
  \item{fixef}{A data.frame with joint effect size estimates (Est), standard errors (SE), chi-squared test statistics (Stat), p-values (pval), and estimated proportion of variance explained (pve) for each of the variants specified in \code{G}.}
  \item{betaCov}{Estimated covariance matrix for the variants in \code{G}.}
}
%\references{
%}

\author{Adrienne M. Stilp, Matthew P. Conomos}
\seealso{
\code{\link{fitNullModel}} for fitting the mixed model and performing the variance component estimation.
\code{\link[GWASTools]{GenotypeData}} and \code{\link[SeqVarTools]{SeqVarData}} for obtaining genotype matrices.
}
\examples{
library(GWASTools)

# File path to GDS file
gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
# Read in GDS data
HapMap_geno <- GdsGenotypeReader(filename = gdsfile)
# Create a GenotypeData class object
HapMap_genoData <- GenotypeData(HapMap_geno)

# Prepare genotypes for the first five SNPs and the first 20 samples.
# Transpose it so that samples are rows and SNPs are columns.
geno <- t(getGenotype(HapMap_genoData, snp = c(1, 5), scan = c(1, 20)))
# Set row and column names.
rownames(geno) <- as.character(GWASTools::getScanID(HapMap_genoData))[1:20]
colnames(geno) <- sprintf("snp\%s", 1:5)

# Create a phenotype
set.seed(5)
phen <- data.frame(
  outcome = rnorm(1:20),
  covar = rnorm(1:20),
  stringsAsFactors = FALSE
)
rownames(phen) <- rownames(geno)

# load saved matrix of KING-robust estimates
data("HapMap_ASW_MXL_KINGmat")
cov.mat = HapMap_ASW_MXL_KINGmat[rownames(phen), rownames(phen)]

# Fit a null model.
nullmod <- fitNullModel(phen,
                        outcome = "outcome",
                        "covars" = "covar",
                        cov.mat = cov.mat)


# Run the joint score test.
jointScoreTest(nullmod, geno)

close(HapMap_genoData)

}
