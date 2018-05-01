\name{nullModelFastSKAT}
\alias{nullModelFastSKAT}

\title{Prepare a null model for fastSKAT}

\description{
  Prepare a null model for fastSKAT
}

\usage{
nullModelFastSKAT(null.model, threshold=1e-10)
}

\arguments{
  \item{null.model}{The output of \code{fitNullModel}.}
  \item{threshold}{Threshold for setting values of the covariance matrix to 0.}
}

\details{
  Computes matrices needed for calling \code{\link{assocTestAggregate}} with \code{test="fastSKAT"}.
}
  
\value{
  An object of class '\code{GENESIS.nullModelFastSKAT}'.
}
  
\references{
  Lumley, T.S. et al. (2018). Sequence kernel association tests for large sets of markers: tail probabilities for large quadratic forms. Submitted.
}

\author{Tamar Sofer, Thomas Lumley}

\seealso{
    \code{\link{fitNullModel}}, \code{\link{assocTestAggregate}}
}
  
\examples{
library(SeqVarTools)
library(Biobase)

# open a sequencing GDS file
gdsfile <- seqExampleFileName("gds")
gds <- seqOpen(gdsfile)

# simulate some phenotype data
data(pedigree)
pedigree <- pedigree[match(seqGetData(gds, "sample.id"), pedigree$sample.id),]
pedigree$outcome <- rnorm(nrow(pedigree))

# construct a SeqVarData object
seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(pedigree))

# fit the null model
nullmod <- fitNullModel(seqData, outcome="outcome", covars="sex")

# prepare the null model for fastSKAT
nullmod <- nullModelFastSKAT(nullmod)

# Window Iterator
seqSetFilterChrom(seqData, include="22")
iterator <- SeqVarWindowIterator(seqData)

# fastSKAT
assoc <- assocTestAggregate(iterator, nullmod, test="fastSKAT")
head(assoc$results)
head(assoc$variantInfo)

seqClose(seqData)
}
\keyword{association}
