\name{match_signif_hits}
\alias{match_signif_hits}

\title{Return Dataframe with Significant Variant Results Matched Across Groups}
\description{
Returns either a list of dataframes or a single dataframe that integrates results of significant variants across groups.
}
\usage{
match_signif_hits(res_list, threshold, return_df=FALSE, variant.id_var=NULL)
}
\arguments{
    \item{res_list}{List of results returned by \code{assocTestSingle_split}.}
    \item{threshold}{Minimum p-value threshold for filtering variants}
    \item{return_df}{If \code{TRUE}, returns a single \code{data.frame} object that combines results across groups. Default is \code{FALSE} and returns a list of \code{data.frame} objects.}
    \item{variant.id_var}{Variable name for variantID in results list if it is not \code{variantID} or \code{variant.id}.}
}

\author{Jiwon Lee}

\seealso{
	\code{\link{find_incomplete_hits}} for creating a vector of variants that were not found in every group after running \code{match_signif_hits}.
     \code{\link{recreate_iterator}} to subset gds file to variants of interest, then format it into an iterator for running single association tests.
     \code{\link{run_split_subset}} for running split carriers only single association tests that return all variants in the iterator, including monomorphic variants.
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

# Create id_list
id_list <- gen_id_list(pedigree, id_var="sample.id", group_var="sex")

# construct a SeqVarIterator object
seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(pedigree))
iterator <- SeqVarBlockIterator(seqData)

# fit the null model
nullmod <- fitNullModel(iterator, outcome="outcome", covars="sex", family="binomial")

# run the association test
assoc <- assocTestSingle_split(iterator, nullmod, id_list=id_list, test="BinomiRare", max.alt.freq=0.5)

# save results for each group separately
save_split_results(assoc, output_prefix='test')

hits <- match_signif_hits(assoc, threshold=0.001)
incomplete <- find_incomplete_hits(hits)

iterator <- recreate_iterator(gds, AnnotatedDataFrame(pedigree), incomplete_variants=incomplete)

sub <- run_split_subset(iterator, nullmod, id_list, test="BinomiRare", keep_all=FALSE)

seqClose(iterator)

}