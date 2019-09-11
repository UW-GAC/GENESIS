\name{mergeNullModelBR}
\alias{mergeNullModelBR}

\title{Combine null models from multiple groups into one}
\description{
Generate a combined null model for use in carriers-only tests.
}
\usage{
mergeNullModelBR(nullmod.list, gdsfile)
}
\arguments{
    \item{nullmod.list}{List of null models to combine}
    \item{gdsfile}{Path to a gds file - this is used to order the sample IDs}
}
\author{Tamar Sofer}

\seealso{
     \code{\link{fitNullModel}} for fitting null models
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

# divide pedigree into two groups
pedigree_one <- pedigree[1:50, ]
pedigree_one_aligned <- data.frame(sample.id=pedigree$sample.id)
pedigree_one_aligned <- dplyr::left_join(pedigree_one_aligned, pedigree_one, by="sample.id")

pedigree_two <- pedigree[51:90, ]
pedigree_two_aligned <- data.frame(sample.id=pedigree$sample.id)
pedigree_two_aligned <- dplyr::left_join(pedigree_two_aligned, pedigree_two, by="sample.id")

# run model for first group
seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(pedigree_one_aligned))
iterator <- SeqVarBlockIterator(seqData)
nullmod_one <- fitNullModel(iterator, outcome="outcome", covars="sex", family="binomial")

# run model for second group
seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(pedigree_two_aligned))
iterator <- SeqVarBlockIterator(seqData)
nullmod_two <- fitNullModel(iterator, outcome="outcome", covars="sex", family="binomial")

# combine null models
nullmod.list <- vector(mode = "list", length = 2)
nullmod.list[[1]] <- nullmod_one
nullmod.list[[2]] <- nullmod_two

seqClose(gds)
gds_path <- system.file("extdata", "CEU_Exon.gds", package="SeqArray")

##combine null models
combined_nullmod <- mergeNullModelBR(nullmod.list,  gds_path)

# setup and  run the association test
gds<- seqExampleFileName("gds")
seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(pedigree))
iterator <- SeqVarBlockIterator(seqData)
assoc <- assocTestSingle(iterator, combined_nullmod, test="BinomiRare", AF.max=0.5)

seqClose(iterator)
}
