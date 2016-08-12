\name{fitNullReg}
\alias{fitNullReg}
\title{Fit a Regression Model Under the Null Hypothesis}
\description{\code{fitNullReg} fits a regression model. The output of \code{fitNullReg} can be passed to \code{\link{assocTestSeq}} or \code{\link{assocTestSeqWindow}} for the purpose of genetic association testing.} 
\usage{
fitNullReg(scanData, outcome, covars = NULL, scan.include = NULL,
            family = gaussian, verbose = TRUE)
}
\arguments{
  \item{scanData}{An object of class \code{ScanAnnotationDataFrame} from the package \code{GWASTools}, \code{AnnotatedDataFrame}, or class \code{data.frame} containing the outcome and covariate data for the samples to be used for the analysis. \code{scanData} must have a column \code{scanID} containing unique IDs for all samples.}
  \item{outcome}{A character string specifying the name of the outcome variable in \code{scanData}.}
  \item{covars}{A vector of character strings specifying the names of the fixed effect covariates in \code{scanData}; an intercept term is automatically included. If NULL (default) the only fixed effect covariate is the intercept term.}
  \item{scan.include}{A vector of scanIDs for samples to include in the analysis.  If NULL, all samples in \code{scanData} are included.}
  \item{family}{A description of the error distribution to be used in the model. The default "gaussian" fits a linear model; see \code{\link{family}} for further options.}
  \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is TRUE.}
}
%\details{
%}
\value{A list including:
    \item{fixef}{A data.frame with effect size estimates (betas), standard errors, chi-squared test statistics, and p-values for each of the fixed effect covariates specified in \code{covars}.}
    \item{betaCov}{The estimated covariance matrix of the effect size estimates (betas) of the fixed effect covariates. This can be used for hypothesis tests regarding the fixed effects.}
    \item{resid.response}{The residuals from the model.}
    \item{logLik}{The log-likelihood value.}
    \item{AIC}{The Akaike Information Criterion value.}
    \item{workingY}{The "working" outcome vector. When \code{family} is gaussian, this is just the original outcome vector. When \code{family} is not gaussian, this is the PQL linearization of the outcome vector. This is used by \code{assocTestSeq} for genetic association testing.}
    \item{model.matrix}{The design matrix for the fixed effect covariates used in the model.}
    \item{aliased}{Coefficients removed from the model.}
    \item{sigma}{Variance of the model.}
    \item{scanID}{A vector of scanIDs for the samples used in the analysis.}
    \item{family}{A character string specifying the family used in the analysis.}
}
%\references{}
\author{Matthew P. Conomos}
%\note{}
%\seealso{
%}
%\examples{
%}
\keyword{association}
