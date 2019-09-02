\name{gen_id_list}
\alias{gen_id_list}

\title{Generate ID list for use in split tests}
\description{
Generate list containing IDs separated by groups for use in \code{assocTestSingle_split} or \code{assocTestAggregate_split}.
}
\usage{
gen_id_list(data, id_var, group_var)
}
\arguments{
    \item{data}{dataframe containing IDs matching genetic data sample IDs, and variable to divide groups by}
    \item{id_var}{name of variable containing IDs}
    \item{group_var}{name of variable to divide groups by}
}
\author{Jiwon Lee}

\seealso{
     \code{\link{assocTestSingle_split}} for performing carriers-only single association tests across different groups
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