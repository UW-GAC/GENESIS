\name{find_incomplete_hits}
\alias{find_incomplete_hits}

\title{Return variants that are not found across all groups}
\description{
Returns vector of variant IDs that are not found across all groups. Meant to be used on list or dataframe of matched results returned by \code{match_signif_hits}, and used for the \code{incomplete_variants} argument in \code{recreate_iterator}.
}
\usage{
find_incomplete_hits(matched_results, n_groups=NULL)
}
\arguments{
    \item{matched_results}{Either a dataframe with a minimum of \code{ref_group} and \code{variantID} columns, or a list containing these dataframes, ideally returned from \code{match_signif_hits}.}
    \item{n_groups}{The minimum number of group results in which variants should be found. Default is \code{NULL}.}
}
\author{Jiwon Lee}

\seealso{
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