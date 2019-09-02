\name{save_split_results}
\alias{save_split_results}

\title{Save split analysis results for each group separately}
\description{
Function to save the list of results returned from \code{assocTestSingle_split} or \code{assocTestAggregate_split} separately by group, given a prefix.
}
\usage{
save_split_results(res_list, output_prefix=NULL)
}
\arguments{
    \item{res_list}{List of results returned from \code{assocTestSingle_split} or \code{assocTestAggregate_split}}
    \item{output_prefix}{Prefix for output file names - file name structure will be prefix_groupname.RData format if provided. Otherwise, it will have a groupname_results.RData format.}
}
\author{Jiwon Lee}

\seealso{
     \code{\link{gen_id_list}} for performing  carriers-only genotype association tests across different groups
  using a single null model
     \code{\link{assocTestAggregate_split}} for performing carriers-only aggregate association tests across different groups using a single null model
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

seqClose(iterator)
}