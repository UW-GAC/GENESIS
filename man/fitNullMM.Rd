\name{fitNullMM}
\alias{fitNullMM}
\title{Fit a Mixed Model Under the Null Hypothesis}
\description{\code{fitNullMM} fits a mixed model with random effects specified by their covariance structures; this allows for the inclusion of a polygenic random effect using a kinship matrix or genetic relationship matrix (GRM). The output of \code{fitNullMM} can be used to estimate genetic heritability and can be passed to \code{\link{assocTestMM}} for the purpose of genetic association testing.} 
\usage{
fitNullMM(scanData, outcome, covars = NULL, covMatList, scan.include = NULL,
            family = gaussian, group.var = NULL, start = NULL, 
            AIREML.tol = 1e-6, maxIter = 100, dropZeros = TRUE, verbose = TRUE)
}
\arguments{
  \item{scanData}{An object of class \code{ScanAnnotationDataFrame} from the package \code{GWASTools}, \code{AnnotatedDataFrame}, or class \code{data.frame} containing the outcome and covariate data for the samples to be used for the analysis. \code{scanData} must have a column \code{scanID} containing unique IDs for all samples.}
  \item{outcome}{A character string specifying the name of the outcome variable in \code{scanData}.}
  \item{covars}{A vector of character strings specifying the names of the fixed effect covariates in \code{scanData}; an intercept term is automatically included. If NULL (default) the only fixed effect covariate is the intercept term.}
  \item{covMatList}{A list of matrices specifying the covariance structures of the random effects terms. The column and row names of each of these matrices must match the scanIDs from \code{scanData}. If only one random effect is being used, a single matrix (not in a list) can be used. See 'Details' for more information.}
  \item{scan.include}{A vector of scanIDs for samples to include in the analysis.  If NULL, all samples in \code{scanData} are included.}
  \item{family}{A description of the error distribution to be used in the model. The default "gaussian" fits a linear mixed model; see \code{\link{family}} for further options, and see 'Details' for more information.}
  \item{group.var}{This variable can only be used when \code{family} = gaussian. A character string specifying the name of a categorical variable in \code{scanData} that is used to fit heterogeneous residual error variances. If NULL (default), then a standard LMM with constant residual variance for all samples is fit. See 'Details' for more information.}
  \item{start}{A vector of starting values for the variance component estimation procedure. The function will pick reasonable starting values when left NULL (default). See 'Details' for more information.}
  \item{AIREML.tol}{The convergence threshold for the Average Information REML (AIREML) procedure used to estimate the variance components of the random effects. See 'Details' for more information.}
  \item{maxIter}{The maximum number of iterations allowed in the AIREML procedure.}
  \item{dropZeros}{Logical indicator of whether variance component terms that converge to 0 should be removed from the model; the default is TRUE. See 'Details' for more information.}
  \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is TRUE.}
}
\details{
    \code{covMatList} is used to specify the covariance structures of the random effects terms in the model.  For example, to include a polygenic random effect, one matrix in \code{covMatList} could be a kinship matrix or a genetic relationship matrix (GRM).  As another example, to include household membership as a random effect, one matrix in \code{covMatList} should be a 0/1 matrix with a 1 in the [i,j] and [j,i] entries if individuals i and j are in the same household and 0 otherwise; the diagonals of such a matrix should all be 1.
    
    When \code{family} is not gaussian, the penalized quasi-likelihood (PQL) approximation to the generalized linear mixed model (GLMM) is fit following the procedure of GMMAT (Chen et al.).
    
    For some outcomes, there may be evidence that different groups of observations have different residual variances, and the standard LMM assumption of homoscedasticity is violated. When \code{group.var} is specified, separate (heterogeneous) residual variance components are fit for each unique value of \code{group.var}.
    
    Let m be the number of matrices in \code{covMatList} and let g be the number of categories in the variable specified by \code{group.var}. The length of the \code{start} vector must be (m + 1) when \code{family} is gaussian and \code{group.var} is NULL; (m + g) when \code{family} is gaussian and \code{group.var} is specified; or m when \code{family} is not gaussian.
    
    A Newton-Raphson iterative procedure with Average Information REML (AIREML) is used to estimate the variance components of the random effects. When the Euclidean distance between the new and previous variance component estimates is less than \code{AIREML.tol}, the algorithm declares convergence of the estimates. Sometimes a variance component may approach the boundary of the parameter space at 0; step-halving is used to prevent any component from becomming negative. However, when a variance component gets near the 0 boundary, the algorithm can sometimes get "stuck", preventing the other variance components from converging; if \code{dropZeros} is TRUE, then variance components that converge to a value less than \code{AIREML.tol} will be dropped from the model and the estimation procedure will continue with the remaining variance components.
}
\value{An object of class '\code{GENESIS.nullMixedModel}'. A list including:
    \item{varComp}{The variance component estimates. There is one variance component for each random effect specified in \code{covMatList}. When \code{family} is gaussian, there are additional residual variance components; one residual variance component when \code{group.var} is NULL, and as many residual variance components as there are unique values of \code{group.var} when it is specified.}
    \item{varCompCov}{The estimated covariance matrix of the variance component estimates given by \code{varComp}. This can be used for hypothesis tests regarding the variance components.}
    \item{fixef}{A data.frame with effect size estimates (betas), standard errors, chi-squared test statistics, and p-values for each of the fixed effect covariates specified in \code{covars}.}
    \item{betaCov}{The estimated covariance matrix of the effect size estimates (betas) of the fixed effect covariates. This can be used for hypothesis tests regarding the fixed effects.}
    \item{fitted.values}{The fitted values from the mixed model; i.e. W*beta where W is the design matrix and beta are the effect size estimates for the fixed effects.}
    \item{resid.marginal}{The marginal residuals from the mixed model; i.e. Y - W*beta where Y is the vector of outcome values.}
    \item{eta}{The linear predictor from the mixed model; i.e. W*beta + Z*u where Z*u specifies the effects of the random effects.}
    \item{resid.conditional}{The conditional residuals from the mixed model; i.e. Y - W*beta - Z*u.}
    \item{logLikR}{The restricted log-likelihood value.}
    \item{logLik}{The log-likelihood value.}
    \item{AIC}{The Akaike Information Criterion value.}
    \item{RSS}{The residual sum of squares from the model fit. When \code{family} is gaussian, this will typically be 1 since the residual variance component is estimated separately.}
    \item{workingY}{The "working" outcome vector. When \code{family} is gaussian, this is just the original outcome vector. When \code{family} is not gaussian, this is the PQL linearization of the outcome vector. This is used by \code{assocTestMM} for genetic association testing. See 'Details' for more information.}
    \item{model.matrix}{The design matrix for the fixed effect covariates used in the model.}
    \item{cholSigmaInv}{The Cholesky decomposition of the inverse of the estimated outcome covariance structure. This is used by \code{assocTestMM} for genetic association testing.}
    \item{scanID}{A vector of scanIDs for the samples used in the analysis.}
    \item{family}{A character string specifying the family used in the analysis.}
    \item{converged}{A logical indicator of whether the AIREML procedure for estimating the random effects variance components converged.}
    \item{zeroFLAG}{A vector of logicals the same length as \code{varComp} specifying whether the corresponding variance component estimate was set to 0 by the function due to convergence to the boundary in the AIREML procedure.}
    \item{hetResid}{A logical indicator of whether heterogeneous residual variance components were used in the model (specified by \code{group.var}).}
}
\references{
    Chen H, Wang C, Conomos MP, Stilp AM, Li Z, Sofer T, Szpiro AA, Chen W, Brehm JM, Celedon JC, Redline S, Papanicolaou GJ, Thornton TA, Laurie CC, Rice K and Lin X. Control for Population Structure and Relatedness for Binary Traits in Genetic Association Studies Using Logistic Mixed Models. (Submitted).
    
    Breslow NE and Clayton DG. (1993). Approximate Inference in Generalized Linear Mixed Models. Journal of the American Statistical Association 88: 9-25.
    
    Gilmour, A.R., Thompson, R., & Cullis, B.R. (1995). Average information REML: an efficient algorithm for variance parameter estimation in linear mixed models. Biometrics, 1440-1450.
    
    Gogarten, S.M., Bhangale, T., Conomos, M.P., Laurie, C.A., McHugh, C.P., Painter, I., ... & Laurie, C.C. (2012). GWASTools: an R/Bioconductor package for quality control and analysis of Genome-Wide Association Studies. Bioinformatics, 28(24), 3329-3331.
}
\author{Matthew P. Conomos}
%\note{}
\seealso{
    \code{\link{varCompCI}} for estimating confidence intervals for the variance components and the proportion of variability (heritability) they explain, \code{\link{assocTestMM}} for running mixed model genetic association tests using the output from \code{fitNullMM}.
    \code{\link[GWASTools:GWASTools-package]{GWASTools}} for a description of the package containing the \code{\link[GWASTools]{ScanAnnotationDataFrame}} class.
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
}
\keyword{mixed model}
\keyword{association}
\keyword{variance component}