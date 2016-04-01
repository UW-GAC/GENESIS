\name{GENESIS-package}
\alias{GENESIS-package}
\alias{GENESIS}
\docType{package}
\title{GENetic EStimation and Inference in Structured samples (GENESIS):  Statistical methods for analyzing genetic data from samples with population structure and/or relatedness}
\description{The GENESIS package provides methodology for estimating, inferring, and accounting for population and pedigree structure in genetic analyses. The current implementation performs PC-AiR (Conomos et al., 2015, Gen Epi) and PC-Relate (Conomos et al., 2016, AJHG). PC-AiR performs a Principal Components Analysis on genome-wide SNP data for the detection of population structure in a sample that may contain known or cryptic relatedness.  Unlike standard PCA, PC-AiR accounts for relatedness in the sample to provide accurate ancestry inference that is not confounded by family structure. PC-Relate uses ancestry representative principal components to adjust for population structure/ancestry and accurately estimate measures of recent genetic relatedness such as kinship coefficients, IBD sharing probabilities, and inbreeding coefficients. Additionally, functions are provided to perform efficient variance component estimation and mixed model association testing for both quantitative and binary phenotypes.}
\details{
\tabular{ll}{
Package: \tab GENESIS\cr
Type: \tab Package\cr
Version: \tab 2.1.7\cr
Date: \tab 2016-4-1\cr
License: \tab GPL-3\cr
Depends: \tab GWASTools\cr
Suggests: \tab gdsfmt, SNPRelate, RUnit, BiocGenerics, knitr\cr
VignetteBuilder: \tab knitr\cr
biocViews: \tab SNP, GeneticVariability, Genetics, StatisticalMethod, DimensionReduction, PrincipalComponent, GenomeWideAssociation, QualityControl, BiocViews\cr

}
The PC-AiR analysis is performed using the \code{\link{pcair}} function, which takes genotype data and pairwise measures of kinship and ancestry divergence as input and returns PC-AiR PCs as the ouput.
The function \code{\link{pcairPartition}} is called within \code{\link{pcair}} and uses the PC-AiR algorithm to partition the sample into an ancestry representative `unrelated subset' and `related subset'.
The function \code{\link{plot.pcair}} can be used to plot pairs of PCs from a class `\code{pcair}' object returned by the function \code{\link{pcair}}.
The function \code{\link{king2mat}} can be used to convert output text files from the KING software (Manichaikul et al., 2010) into an R matrix of pairwise kinship coefficient estimates in a format that can be used by the functions \code{\link{pcair}} and \code{\link{pcairPartition}}.
The PC-Relate analysis is performed using the \code{\link{pcrelate}} function, which takes genotype data and PCs from PC-AiR and returns estimates of kinship coefficients, IBD sharing probabilities, and inbreeding coefficients. The functions \code{\link{pcrelateReadKinship}}, \code{\link{pcrelateReadInbreed}}, and \code{\link{pcrelateMakeGRM}} provide utilities for reading and making tables or matrices of the PC-Relate output.
There are two functions required to perform SNP genotype association testing with mixed models. First, \code{\link{fitNullMM}} is called to fit the null model (i.e. no SNP genotype term) including fixed effects covariates, such as PC-AiR PCs, and random effects specified by their covariance structures, such as a kinship matrix created from PC-Relate output using \code{\link{pcrelateMakeGRM}}. The function \code{\link{fitNullMM}} uses AIREML to estimate variance components for the random effects, and the function \code{\link{varCompCI}} can be used to find confidence intervals on the estimates as well as the proportion of total variability they explain; this allows for heritability estimation. Second, \code{\link{assocTestMM}} is called with the null model output and the genotype data to perform either Wald or score based association tests.
}
\author{
Matthew P. Conomos and Timothy Thornton

Maintainer: Matthew P. Conomos <mconomos@uw.edu>
}
\references{
  Conomos M.P., Reiner A.P., Weir B.S., & Thornton T.A. (2016). Model-free Estimation of Recent Genetic Relatedness. American Journal of Human Genetics, 98(1), 127-148.
  
	Conomos M.P., Miller M., & Thornton T. (2015). Robust Inference of Population Structure for Ancestry Prediction and Correction of Stratification in the Presence of Relatedness. Genetic Epidemiology, 39(4), 276-293.
    
    Gogarten, S. M., Bhangale, T., Conomos, M. P., Laurie, C. A., McHugh, C. P., Painter, I., ... & Laurie, C. C. (2012). GWASTools: an R/Bioconductor package for quality control and analysis of Genome-Wide Association Studies. Bioinformatics, 28(24), 3329-3331.
    
	Manichaikul, A., Mychaleckyj, J.C., Rich, S.S., Daly, K., Sale, M., & Chen, W.M. (2010). Robust relationship inference in genome-wide association studies. Bioinformatics, 26(22), 2867-2873.

}
%\keyword{ package }
%\seealso{
% ~~ Optional links to other man pages, e.g. ~~
% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
% ~~ simple examples of the most important functions ~~
%}

