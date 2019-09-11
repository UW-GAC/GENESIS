\name{runSplitSubset}
\alias{runSplitSubset}


\title{Run Split Carriers Only Tests that Return All Results Across Groups}

\description{\code{runSplitSubset} performs carriers-only genotype association tests across different groups
  using a single null model fit with \code{\link{fitNullModel}}. It is meant to be run on a subset of variants for which
  results were not returned in at least one group. Monomorphic variants will also be returned with \code{NA}.}

\usage{
runSplitSubset(gdsobj, null.model, id.list, test=c("BinomiRare", "CMP"), remove.groups=NULL,
                   sparse=TRUE, imputed=FALSE, keep.all=TRUE, 
                   male.diploid=TRUE,  genome.build=c("hg19", "hg38"), verbose=TRUE) 
}

\arguments{
    \item{gdsobj}{An object of class \code{\link{SeqVarIterator}} from the package \pkg{\link{SeqVarTools}}, or an object of class \code{\link{GenotypeIterator}} from the package \pkg{\link{GWASTools}}, containing the genotype data for the variants and samples to be used for the analysis.}
    \item{null.model}{A null model object returned by \code{\link{fitNullModel}}.}
    \item{id.list}{List that includes vectors of sample IDs divided by group}
    \item{test}{A character string specifying the type of carriers-only test to be performed. Options are \code{"BinomiRare"} and \code{"CMP"}.}
    \item{remove.groups}{Vector of groups to remove from id.list}
    \item{sparse}{Logical indicator of whether to read genotypes as sparse Matrix objects; the default is \code{TRUE}. Set this to \code{FALSE} if the alternate allele dosage of the genotypes in the test are not expected to be mostly 0.}
    \item{imputed}{Logical indicator of whether to read dosages from the "DS" field containing imputed dosages instead of counting the number of alternate alleles.}
    \item{keep.all}{Logical indicator of whether the overall sample should be tested as well, in addition to each group; default is \code{TRUE}.}
    \item{male.diploid}{Logical for whether males on sex chromosomes are coded as diploid.}
    \item{genome.build}{A character sting indicating genome build; used to identify pseudoautosomal regions on the X and Y chromosomes.}
    \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is \code{TRUE}.}

}

\details{
    All variants including will be returned, for conveniently comparing results across groups. 
    To simplify finding and rerunning tests on variants of interest, reference \code{matchSignifHits}, 
    \code{findIncompleteHits}, and \code{recreateIterator}. 
    Group-specific results can be saved as a separate RData file by calling \code{saveSplitResults}. This function is currently only compatible with \code{SeqVarIterator} objects.
}

\value{A list of results by group. Each set of results is a data.frame where each row refers to a different variant with the columns:
    \item{variant.id}{The variant ID}
    \item{chr}{The chromosome value}
    \item{pos}{The base pair position}
    \item{allele.index}{The index of the alternate allele. For biallelic variants, this will always be 1.}
    \item{n.obs}{The number of samples with non-missing genotypes}
    \item{freq}{The estimated alternate allele frequency}
    If \code{test} is \code{"BinomiRare"}:
    \item{n.carrier}{Number of individuals with at least one copy of the alternate allele}
    \item{n.D.carrier}{Number of cases with at least one copy of the alternate allele}
    \item{pval}{BinomiRare p-value}
    If \code{test} is \code{"CMP"}:
    \item{n.carrier}{Number of individuals with at least one copy of the alternate allele}
    \item{n.D.carrier}{Number of cases with at least one copy of the alternate allele}
    \item{pval}{CMP p-value}
    \item{mid.pval}{CMP mid-p-value}
}

\author{Jiwon Lee, Tamar Sofer, Matthew P. Conomos, Stephanie M. Gogarten, Ken Rice, Chaoyu Yu}

\seealso{
     \code{\link{matchSignifHits}} for creating a list or dataframe of significant variants that integrates results across groups.
    \code{\link{findIncompleteHits}} for creating a vector of variants that were not found in every group after running \code{matchSignifHits}.
    \code{\link{recreateIterator}} for recreating the iterator object that subsets the gds file to only a subset of variants (called by the variant IDs). Meant to be used in conjunction with \code{findIncompleteHits}.
     \code{\link{saveSplitResults} for saving results list as separate files by group}
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
\keyword{mixed model}
\keyword{association}
