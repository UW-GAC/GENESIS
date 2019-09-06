\name{recreate_iterator}
\alias{recreate_iterator}

\title{Recreate Iterator with Subset of Variants to Run}
\description{
Subset gds file to variants of interest, then format it into an iterator for running single association tests. This was specifically written for use with \code{run_split_subset}, but can be used as a convenience function for other single association tests as well.
}
\usage{
recreate_iterator(gds, annot, incomplete_variants, block.size=1024)
}
\arguments{
    \item{gds}{gds file connection called by \code{seqOpen}}
    \item{annot}{Annotated data frame}
    \item{incomplete_variants}{Vector of variant IDs that should be included in the iterator}
    \item{block.size}{Number of variants to include in each iteration}
}
\author{Jiwon Lee}

\seealso{
     \code{\link{find_incomplete_hits}} for finding variants from split results that are not found across all groups. The returned vector can be directly used for the \code{incomplete_variants} argument.
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