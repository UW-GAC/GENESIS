\name{admixMapMM}
\alias{admixMapMM}

\title{admixMapMM}
\description{
Run admixture analyses
}
\usage{
admixMapMM(admixDataList, nullMMobj, snp.include = NULL,
           chromosome = NULL, snp.block.size = 5000,
           verbose = TRUE)
}

\arguments{
  \item{admixDataList}{named list of \code{\link{GenotypeData}} objects}
    \item{nullMMobj}{A null model object returned by \code{fitNullMM}.}
    \item{snp.include}{A vector of SNP IDs to include in the analysis.  If NULL, see \code{chromosome} for further details.}
    \item{chromosome}{A vector of integers specifying which chromosomes to analyze.  This parameter is only considered when \code{snp.include} is NULL; if \code{chromosome} is also NULL, then all SNPs are included.}
    \item{snp.block.size}{The number of SNPs to read-in/analyze at once.  The default value is 5000.}
    \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is TRUE.}
}


\value{data frame with admixture mapping results}

\author{Matt Conomos}

\seealso{\code{\link{GenotypeData}}, \code{\link{fitNullMM}}, \code{\link{assocTestMM}}}

