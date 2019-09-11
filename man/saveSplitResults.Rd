\name{saveSplitResults}
\alias{saveSplitResults}

\title{Save split analysis results for each group separately}
\description{
Function to save the list of results returned from \code{assocTestSingleSplit} or \code{assocTestAggregateSplit} separately by group, given a prefix.
}
\usage{
saveSplitResults(res.list, output.prefix=NULL)
}
\arguments{
    \item{res.list}{List of results returned from \code{assocTestSingleSplit} or \code{assocTestAggregateSplit}}
    \item{output.prefix}{Prefix for output file names - file name structure will be prefix_groupname.RData format if provided. Otherwise, it will have a groupname_results.RData format.}
}
\author{Jiwon Lee}

\seealso{
     \code{\link{genIDList}} for performing  carriers-only genotype association tests across different groups
  using a single null model
     \code{\link{assocTestAggregateSplit}} for performing carriers-only aggregate association tests across different groups using a single null model
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
pedigree$outcome <- rbinom(nrow(pedigree), 1, 0.4)

# Create id.list
id.list <- genIDList(pedigree, id.var="sample.id", group.var="sex")

# construct a SeqVarIterator object
seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(pedigree))
iterator <- SeqVarBlockIterator(seqData)

# fit the null model
nullmod <- fitNullModel(iterator, outcome="outcome", covars="sex", family="binomial")

# run the association test
assoc <- assocTestSingleSplit(iterator, nullmod, id.list=id.list, test="BinomiRare", AF.max=0.5)

# save results for each group separately
saveSplitResults(assoc, output.prefix='test')

seqClose(iterator)
}