\name{findIncompleteHits}
\alias{findIncompleteHits}

\title{Return variants that are not found across all groups}
\description{
Returns vector of variant IDs that are not found across all groups. Meant to be used on list or dataframe of matched results returned by \code{matchSignifHits}, and used for the \code{incomplete.variants} argument in \code{recreateIterator}.
}
\usage{
findIncompleteHits(matched.results, n.groups=NULL, variant.id.var=NULL)
}
\arguments{
    \item{matched.results}{Either a dataframe with a minimum of \code{ref_group} and \code{variantID} columns, or a list containing these dataframes, ideally returned from \code{matchSignifHits}.}
    \item{n.groups}{The minimum number of group results in which variants should be found. Default is \code{NULL}.}
    \item{variant.id.var}{Variable name for variantID in results list if it is not \code{variantID} or \code{variant.id}.}
}
\author{Jiwon Lee}

\seealso{
     \code{\link{recreateIterator}} to subset gds file to variants of interest, then format it into an iterator for running single association tests.
     \code{\link{runSplitSubset}} for running split carriers only single association tests that return all variants in the iterator, including monomorphic variants.
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

hits <- matchSignifHits(assoc, threshold=0.1)
incomplete <- findIncompleteHits(hits)

iterator <- recreateIterator(gds, AnnotatedDataFrame(pedigree), incomplete.variants=incomplete)

sub <- runSplitSubset(iterator, nullmod, id.list, test="BinomiRare", keep.all=FALSE)

seqClose(iterator)

}