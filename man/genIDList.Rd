\name{genIDList}
\alias{genIDList}

\title{Generate ID list for use in split tests}
\description{
Generate list containing IDs separated by groups for use in \code{assocTestSingleSplit} or \code{assocTestAggregateSplit}.
}
\usage{
genIDList(data, id.var, group.var)
}
\arguments{
    \item{data}{dataframe containing IDs matching genetic data sample IDs, and variable to divide groups by}
    \item{id.var}{name of variable containing IDs}
    \item{group.var}{name of variable to divide groups by}
}
\author{Jiwon Lee}

\seealso{
     \code{\link{assocTestSingleSplit}} for performing carriers-only single association tests across different groups
  using a single null model. 
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