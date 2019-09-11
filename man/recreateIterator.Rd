\name{recreateIterator}
\alias{recreateIterator}

\title{Recreate Iterator with Subset of Variants to Run}
\description{
Subset gds file to variants of interest, then format it into an iterator for running single association tests. This was specifically written for use with \code{runSplitSubset}, but can be used as a convenience function for other single association tests as well.
}
\usage{
recreateIterator(gds, annot, incomplete.variants, block.size=1024)
}
\arguments{
    \item{gds}{gds file connection called by \code{seqOpen}}
    \item{annot}{Annotated data frame}
    \item{incomplete.variants}{Vector of variant IDs that should be included in the iterator}
    \item{block.size}{Number of variants to include in each iteration}
}
\author{Jiwon Lee}

\seealso{
     \code{\link{findIncompleteHits}} for finding variants from split results that are not found across all groups. The returned vector can be directly used for the \code{incomplete.variants} argument.
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
saveSplitResults(assoc, output_prefix='test')

hits <- matchSignifHits(assoc, threshold=0.001)
incomplete <- findIncompleteHits(hits)

iterator <- recreateIterator(gds, AnnotatedDataFrame(pedigree), incomplete.variants=incomplete)

sub <- runSplitSubset(iterator, nullmod, id.list, test="BinomiRare", keep.all=FALSE)

seqClose(iterator)

}