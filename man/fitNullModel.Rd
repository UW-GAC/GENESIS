\name{fitNullModel}
\alias{fitNullModel}
\alias{fitNullModel-methods}
\alias{fitNullModel,AnnotatedDataFrame-method}
\alias{fitNullModel,data.frame-method}
\alias{fitNullModel,SeqVarData-method}
\alias{nullModelInvNorm}

\title{Fit a Model Under the Null Hypothesis}

\description{
  \code{fitNullModel} fits a regression model or a mixed
  model with random effects specified by their covariance structures;
  this allows for the inclusion of a polygenic random effect using a
  kinship matrix or genetic relationship matrix (GRM). The output of
  \code{fitNullModel} can be used to estimate genetic heritability and
  can be passed to \code{\link{assocTestSingle}} or
  \code{\link{assocTestAggregate}} for the purpose of genetic
  association testing.

  \code{nullModelInvNorm} does an inverse normal transform of a previously fit null model.
}

\usage{
\S4method{fitNullModel}{data.frame}(x, outcome, covars = NULL, cov.mat = NULL,
            group.var = NULL, family = "gaussian", start = NULL, 
            AIREML.tol = 1e-6, max.iter = 100, drop.zeros = TRUE, verbose = TRUE)
\S4method{fitNullModel}{AnnotatedDataFrame}(x, outcome, covars = NULL, cov.mat = NULL,
            group.var = NULL, sample.id = NULL, ...)
\S4method{fitNullModel}{SeqVarData}(x, ...)

nullModelInvNorm(null.model, cov.mat = NULL, norm.option = c("by.group", "all"),
                 rescale = c("none", "model", "residSD"),
                 AIREML.tol = 1e-6, max.iter = 100, verbose = TRUE)
}

\arguments{
  \item{x}{An object of class \code{data.frame}, \code{\link{AnnotatedDataFrame}}, or \code{\link{SeqVarData}} containing the outcome and covariate data for the samples to be used for the analysis.}
  \item{outcome}{A character string specifying the name of the outcome variable in \code{x}.}
  \item{covars}{A vector of character strings specifying the names of the fixed effect covariates in \code{x}; an intercept term is automatically included. If \code{NULL} (default) the only fixed effect covariate is the intercept term.}
  \item{cov.mat}{A matrix or list of matrices specifying the covariance structures of the random effects terms. Objects from the \pkg{\link{Matrix}} package are supported. See 'Details' for more information.}
  \item{group.var}{This variable can only be used when \code{family = "gaussian"}. A character string specifying the name of a categorical variable in \code{x} that is used to fit heterogeneous residual error variances. If \code{NULL} (default), then a standard LMM with constant residual variance for all samples is fit. See 'Details' for more information.}
  \item{sample.id}{A vector of IDs for samples to include in the analysis.  If \code{NULL}, all samples in \code{x} are included.}
  \item{family}{A description of the error distribution to be used in the model. The default \code{"gaussian"} fits a linear model; see \code{\link{family}} for further options, and see 'Details' for more information.}
  \item{start}{A vector of starting values for the variance component estimation procedure. The function will pick reasonable starting values when left \code{NULL} (default). See 'Details' for more information.}
  \item{AIREML.tol}{The convergence threshold for the Average Information REML (AIREML) procedure used to estimate the variance components of the random effects. See 'Details' for more information.}
  \item{max.iter}{The maximum number of iterations allowed in the AIREML procedure.}
  \item{drop.zeros}{Logical indicator of whether variance component terms that converge to 0 should be removed from the model; the default is TRUE. See 'Details' for more information.}
  \item{verbose}{Logical indicator of whether updates from the function
    should be printed to the console; the default is TRUE.}
  \item{...}{Arguments to pass to other methods.}
  \item{null.model}{The output of \code{fitNullModel}.}
  \item{norm.option}{Whether the normalization should be done separately within each value of \code{group.var} (\code{"by.group"}) or with all samples together (\code{"all"}).}
  \item{rescale}{Controls whether to rescale the variance for each group after inverse-normal transform, restoring it to the original variance before the transform. \code{"none"} for no rescaling of the residuals; \code{"model"} for model-based rescaling, and \code{"residSD"} to rescale to the standard deviation of the marginal residuals.}
}

\details{
    \code{cov.mat} is used to specify the covariance structures of the random effects terms in the model.  For example, to include a polygenic random effect, one matrix in \code{cov.mat} could be a kinship matrix or a genetic relationship matrix (GRM).  As another example, to include household membership as a random effect, one matrix in \code{cov.mat} should be a 0/1 matrix with a 1 in the \code{[i,j]} and \code{[j,i]} entries if individuals \code{i} and \code{j} are in the same household and 0 otherwise; the diagonals of such a matrix should all be 1.
    
    When \code{family} is not gaussian, the penalized quasi-likelihood (PQL) approximation to the generalized linear mixed model (GLMM) is fit following the procedure of GMMAT (Chen et al.).
    
    For some outcomes, there may be evidence that different groups of observations have different residual variances, and the standard LMM assumption of homoscedasticity is violated. When \code{group.var} is specified, separate (heterogeneous) residual variance components are fit for each unique value of \code{group.var}.
    
    Let \code{m} be the number of matrices in \code{cov.mat} and let \code{g} be the number of categories in the variable specified by \code{group.var}. The length of the \code{start} vector must be \code{(m + 1)} when \code{family} is gaussian and \code{group.var} is \code{NULL}; \code{(m + g)} when \code{family} is gaussian and \code{group.var} is specified; or m when \code{family} is not gaussian.
    
    A Newton-Raphson iterative procedure with Average Information REML (AIREML) is used to estimate the variance components of the random effects. When the Euclidean distance between the new and previous variance component estimates is less than \code{AIREML.tol}, the algorithm declares convergence of the estimates. Sometimes a variance component may approach the boundary of the parameter space at 0; step-halving is used to prevent any component from becomming negative. However, when a variance component gets near the 0 boundary, the algorithm can sometimes get "stuck", preventing the other variance components from converging; if \code{drop.zeros} is TRUE, then variance components that converge to a value less than \code{AIREML.tol} will be dropped from the model and the estimation procedure will continue with the remaining variance components.
}
  
\value{An object of class '\code{GENESIS.nullModel}' or '\code{GENESIS.nullMixedModel}'. A list including:
    \item{family}{A character string specifying the family used in the analysis.}
    \item{hetResid}{A logical indicator of whether heterogeneous residual variance components were used in the model (specified by \code{group.var}).}
    \item{varComp}{The variance component estimates. There is one variance component for each random effect specified in \code{cov.mat}. When \code{family} is gaussian, there are additional residual variance components; one residual variance component when \code{group.var} is \code{NULL}, and as many residual variance components as there are unique values of \code{group.var} when it is specified.}
    \item{varCompCov}{The estimated covariance matrix of the variance component estimates given by \code{varComp}. This can be used for hypothesis tests regarding the variance components.}
    \item{fixef}{A data.frame with effect size estimates (betas), standard errors, chi-squared test statistics, and p-values for each of the fixed effect covariates specified in \code{covars}.}
    \item{betaCov}{The estimated covariance matrix of the effect size estimates (betas) of the fixed effect covariates. This can be used for hypothesis tests regarding the fixed effects.}
    \item{fitted.values}{The fitted values from the model; i.e. W*beta where W is the design matrix and beta are the effect size estimates for the fixed effects.}
    \item{resid.marginal}{The marginal residuals from the model; i.e. Y - W*beta where Y is the vector of outcome values.}
    \item{resid.conditional}{The conditional residuals from the model; i.e. Y - W*beta - Z*u.}
    \item{logLik}{The log-likelihood value.}
    \item{logLikR}{The restricted log-likelihood value.}
    \item{AIC}{The Akaike Information Criterion value.}
    \item{workingY}{The "working" outcome vector. When \code{family} is gaussian, this is just the original outcome vector. When \code{family} is not gaussian, this is the PQL linearization of the outcome vector. This is used by \code{\link{assocTestSingle}} or \code{\link{assocTestAggregate}} for genetic association testing. See 'Details' for more information.}
    \item{outcome}{The original outcome vector, as a 1-column matrix with column name. When \code{family} is gaussian, this is equal to \code{workingY}.}
    \item{model.matrix}{The design matrix for the fixed effect
      covariates used in the model.}
    \item{group.idx}{If \code{group.var} is not \code{NULL}, a list of
      indices for samples in each group.}
    \item{cholSigmaInv}{The Cholesky decomposition of the inverse of the estimated outcome covariance structure. This is used by \code{\link{assocTestSingle}} or \code{\link{assocTestAggregate}} for genetic association testing.}
    \item{converged}{A logical indicator of whether the AIREML procedure for estimating the random effects variance components converged.}
    \item{zeroFLAG}{A vector of logicals the same length as \code{varComp} specifying whether the corresponding variance component estimate was set to 0 by the function due to convergence to the boundary in the AIREML procedure.}
    \item{RSS}{The residual sum of squares from the model fit. When \code{family} is gaussian, this will typically be 1 since the residual variance component is estimated separately.}
    \item{sample.id}{A vector of IDs for the samples used in the analysis.}
}
  
\references{
    Chen H, Wang C, Conomos MP, Stilp AM, Li Z, Sofer T, Szpiro AA, Chen
    W, Brehm JM, Celedon JC, Redline S, Papanicolaou GJ, Thornton TA,
    Laurie CC, Rice K and Lin X. (2016) Control for Population Structure and
    Relatedness for Binary Traits in Genetic Association Studies Using
    Logistic Mixed Models. American Journal of Human Genetics, 98(4):653-66.
    
    Breslow NE and Clayton DG. (1993). Approximate Inference in Generalized Linear Mixed Models. Journal of the American Statistical Association 88: 9-25.
    
    Gilmour, A.R., Thompson, R., & Cullis, B.R. (1995). Average information REML: an efficient algorithm for variance parameter estimation in linear mixed models. Biometrics, 1440-1450.
}

\author{Matthew P. Conomos, Stephanie M. Gogarten, Tamar Sofer, Ken Rice, Chaoyu Yu}

\seealso{
    \code{\link{varCompCI}} for estimating confidence intervals for the variance components and the proportion of variability (heritability) they explain, \code{\link{assocTestSingle}} or \code{\link{assocTestAggregate}} for running genetic association tests using the output from \code{fitNullModel}.
}
  
\examples{
library(GWASTools)
library(GENESIS)

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

annot <- data.frame(sample.id = mypcrel$sample.id, 
                    pc1 = mypcair$vectors[,1], pheno = pheno)

# make cov.mat
cov.mat <- pcrelateMakeGRM(mypcrel)

# fit the null mixed model
nullmod <- fitNullModel(annot, outcome = "pheno", covars = "pc1", cov.mat = cov.mat)
}
\keyword{mixed model}
\keyword{association}
\keyword{variance component}
