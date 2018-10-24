\name{kin2gds}
\alias{kin2gds}
\alias{mat2gds}
\title{Store kinship matrix in GDS}
\description{\code{kin2gds} and \code{mat2gds} write kinship matrices to
  GDS files.
}
\usage{
kin2gds(ibdobj, gdsfile)
mat2gds(mat, gdsfile)
}
\arguments{
  \item{ibdobj}{A list with elements \code{sample.id} and \code{kinship}, such as the output of \code{\link[SNPRelate]{snpgdsIBDKING}}.}
  \item{mat}{A matrix with sample identifiers in colnames.}
  \item{gdsfile}{A character string with the name of the GDS file to create.}
}
\details{
  When using \code{pcair} or \code{pcairPartition} with large sample sizes, it can be more memory efficient to read data from GDS files. \code{kin2gds} and \code{mat2gds} store kinship matrices in GDS files for use with these functions.
}
\author{Stephanie M. Gogarten}
\seealso{\code{\link{kingToMatrix}}, \code{\link[SNPRelate]{snpgdsIBDKING}}}
\examples{
library(gdsfmt)

# KING results from the command-line program
file.kin0 <- system.file("extdata", "MXL_ASW.kin0", package="GENESIS")
file.kin <- system.file("extdata", "MXL_ASW.kin", package="GENESIS")
KINGmat <- kingToMatrix(c(file.kin0, file.kin))
gdsfile <- tempfile()
mat2gds(KINGmat, gdsfile)
gds <- openfn.gds(gdsfile)
gds
closefn.gds(gds)

# KING results from SNPRelate
library(SNPRelate)
geno <- snpgdsOpen(snpgdsExampleFileName())
king <- snpgdsIBDKING(geno)
closefn.gds(geno)
kin2gds(king, gdsfile)
gds <- openfn.gds(gdsfile)
gds
closefn.gds(gds)
}