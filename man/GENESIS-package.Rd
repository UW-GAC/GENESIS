\name{GENESIS-package}
\alias{GENESIS-package}
\alias{GENESIS}
\docType{package}
\title{GENetic EStimation and Inference in Structured samples (GENESIS):  Statistical methods for analyzing genetic data from samples with population structure and/or relatedness}
\description{The GENESIS package provides methodology for estimating, inferring, and accounting for population and pedigree structure in genetic analyses.  The current implementation performs PC-AiR (Conomos et al., 2015): a Principal Components Analysis on genome-wide SNP data for the detection of population structure in a sample that may contain known or cryptic relatedness.  Unlike standard PCA, PC-AiR accounts for relatedness in the sample to provide accurate ancestry inference that is not confounded by family structure.}
\details{
\tabular{ll}{
Package: \tab GENESIS\cr
Type: \tab Package\cr
Version: \tab 0.99.4\cr
Date: \tab 2015-03-23\cr
License: \tab GPL-3\cr
Depends: \tab GWASTools\cr
Suggests: \tab gdsfmt, SNPRelate, RUnit, BiocGenerics, knitr\cr
VignetteBuilder: \tab knitr\cr
biocViews: \tab SNP, GeneticVariability, Genetics, StatisticalMethod, DimensionReduction, PrincipalComponent, GenomeWideAssociation, QualityControl, BiocViews\cr

}
The main function in this package is \code{\link{pcair}}.  This function takes genotype data and pairwise measures of kinship and ancestry divergence as input and returns PC-AiR PCs as the ouput.
The function \code{\link{pcairPartition}} is called within \code{\link{pcair}} and uses the PC-AiR algorithm to partition the sample into an ancestry representative `unrelated subset' and `related subset'.
The function \code{\link{plot.pcair}} can be used to plot pairs of PCs from a class `\code{pcair}' object returned by the function \code{\link{pcair}}.
The function \code{\link{king2mat}} can be used to convert output text files from the KING software (Manichaikul et al., 2010) into an R matrix of pairwise kinship coefficient estimates in a format that can be used by the functions \code{\link{pcair}} and \code{\link{pcairPartition}}.
}
\author{
Matthew P. Conomos and Timothy Thornton

Maintainer: Matthew P. Conomos <mconomos@uw.edu>
}
\references{
	Conomos M.P., Miller M., & Thornton T. Robust Inference of Population Structure for Ancestry Prediction and Correction of Stratification in the Presence of Relatedness. (Accepted to Genetic Epidemiology).
    
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

