\name{effectAllele}
\alias{effectAllele}
\alias{effectAllele-methods}
\alias{effectAllele,SeqVarGDSClass-method}
\alias{effectAllele,GenotypeData-method}
\title{Return the effect allele for association testing}
\description{\code{effectAllele} returns the effect allele for association testing.}
\usage{
\S4method{effectAllele}{SeqVarGDSClass}(gdsobj, variant.id=NULL)
\S4method{effectAllele}{GenotypeData}(gdsobj, variant.id=NULL)
}
\arguments{
  \item{gdsobj}{An object of class \code{\link{SeqVarIterator}} from the package \pkg{\link{SeqVarTools}}, or an object of class \code{\link{GenotypeIterator}} from the package \pkg{\link{GWASTools}}, containing the genotype data for the variants and samples to be used for the analysis.}
  \item{variant.id}{A vector of identifiers for variants to return.}
}
\details{
  \code{effectAllele} returns the effect allele corresponding to association test results from \code{\link{assocTestSingle}} or \code{\link{assocTestAggregate}}. \code{variant.id} allows the user to specify for which variants effect alleles should be returned.
}
\value{
  A data.frame with the following columns:
  \item{variant.id}{The variant ID}
  \item{effect.allele}{The character value for the effect allele}
  \item{other.allele}{The character value for the other (non-effect) allele}
}

\author{Stephanie M. Gogarten}
\seealso{
\code{\link{assocTestSingle}}, \code{\link{assocTestAggregate}}
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

# add effect allele to the results
eff <- effectAllele(seqData, variant.id=assoc$variant.id)
assoc <- merge(assoc, eff)
head(assoc)

seqClose(iterator)
}
\keyword{mixed model}
\keyword{association}
