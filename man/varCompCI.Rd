\name{varCompCI}
\alias{varCompCI}
\title{Variance Component Confidence Intervals}
\description{\code{varCompCI} provides confidence intervals for the variance component estimates found using \code{\link{fitNullMM}}. The confidence intervals can be found on either the original scale or for the proportion of total variability explained.}
\usage{
varCompCI(nullMMobj, prop = TRUE)
}
\arguments{
  \item{nullMMobj}{A null model object returned by \code{\link{fitNullMM}}.}
  \item{prop}{A logical indicator of whether the point estimates and confidence intervals should be returned as the proportion of total variability explained (TRUE) or on the orginal scale (FALSE).}
}
\details{
    \code{varCompCI} takes the object returned by \code{\link{fitNullMM}} as its input and returns point estimates and confidence intervals for each of the random effects variance component estimates. If a kinship matrix or genetic relationship matrix (GRM) was included as a random effect in the model fit using \code{\link{fitNullMM}}, then this function can be used to provide a heritability estimate when \code{prop} is TRUE.
}
\value{\code{varCompCI} prints a table of point estimates and 95\% confidence interval limits for each estimated variance component.}
%\references{
%}
\author{Matthew P. Conomos}
\seealso{\code{\link{fitNullMM}} for fitting the mixed model and performing the variance component estimation.}
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

# find the variance component CIs
varCompCI(nullmod, prop = TRUE)
varCompCI(nullmod, prop = FALSE)
}
\keyword{variance component}
\keyword{mixed model}
\keyword{heritability}
