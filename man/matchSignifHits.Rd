\name{matchSignifHits}
\alias{matchSignifHits}

\title{Return Dataframe with Significant Variant Results Matched Across Groups}
\description{
Returns either a list of dataframes or a single dataframe that integrates results of significant variants across groups.
}
\usage{
matchSignifHits(res.list, threshold, return.df=FALSE, variant.id.var=NULL)
}
\arguments{
    \item{res.list}{List of results returned by \code{assocTestSingleSplit}.}
    \item{threshold}{Minimum p-value threshold for filtering variants}
    \item{return.df}{If \code{TRUE}, returns a single \code{data.frame} object that combines results across groups. Default is \code{FALSE} and returns a list of \code{data.frame} objects.}
    \item{variant.id.var}{Variable name for variantID in results list if it is not \code{variantID} or \code{variant.id}.}
}

\author{Jiwon Lee}

\seealso{
	\code{\link{findIncompleteHits}} for creating a vector of variants that were not found in every group after running \code{matchSignifHits}.
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