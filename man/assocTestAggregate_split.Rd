\name{assocTestAggregate_split}
\alias{assocTestAggregate_split}
\alias{assocTestAggregate_split-methods}
\alias{assocTestAggregate_split,SeqVarIterator-method}

\title{Aggregate Association Testing}
\description{\code{assocTestAggregate_split} performs carriers-only aggregate association tests across different groups using a single null model fit with \code{\link{fitNullModel}}.}

\usage{
\S4method{assocTestAggregate_split}{SeqVarIterator}(gdsobj, null.model, id_list, AF.max=1,
                   burden.test=c("BinomiRare", "CMP"), keep_all=TRUE,
                   sparse=TRUE, imputed=FALSE, male.diploid=TRUE, genome.build=c("hg19", "hg38"), verbose=TRUE)
}

\arguments{
    \item{gdsobj}{An object of class \code{\link{SeqVarIterator}} from the package \pkg{\link{SeqVarTools}} containing the genotype data for the variants and samples to be used for the analysis.}
    \item{null.model}{A null model object returned by \code{fitNullModel}.}
    \item{id_list}{List that includes vectors of sample IDs divided by group}
    \item{AF.max}{A numeric value specifying the upper bound on the alternate allele frequency for variants to be included in the analysis. The default is \code{1} (include all), but it is recommended to set a low threshold.}
    \item{burden.test}{A character string specifying which type of carriers-onyl test to perform. Options are \code{"BinomiRare"}, and \code{"CMP"}.}
     \item{keep_all}{Logical indicator of whether the overall sample should be tested as well, in addition to each group; default is \code{TRUE}.}
    \item{sparse}{Logical indicator of whether to read genotypes as sparse Matrix objects; the default is \code{TRUE}. Set this to \code{FALSE} if the alternate allele dosage of the genotypes in the test are not expected to be mostly 0.}
    \item{imputed}{Logical indicator of whether to read dosages from the "DS" field containing imputed dosages instead of counting the number of alternate alleles.}
    \item{male.diploid}{Logical for whether males on sex chromosomes are coded as diploid.}
    \item{genome.build}{A character sting indicating genome build; used to identify pseudoautosomal regions on the X and Y chromosomes.}
    \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is \code{TRUE}.}
}

\details{
    The type of aggregate unit tested depends on the class of iterator used for \code{gdsobj}. Options include sliding windows, specific ranges of variants or selection of individual variants (ranges with width 1). See \code{\link{SeqVarIterator}} for more details.
  
    The effect size estimate is for each copy of the alternate allele.
    For multiallelic variants, each alternate allele is tested separately.

    This function is currently only compatible with SeqVarIterator objects.
}

\value{A list with the following items:
  \item{results}{A data.frame containing the results from the main analysis. Each row is a separate aggregate test:}
    If \code{gdsobj} is a \code{\link{SeqVarWindowIterator}}:
    \item{chr}{The chromosome value}
    \item{start}{The start position of the window}
    \item{end}{The end position of the window}

    Always:
    \item{n.site}{The number of variant sites included in the test.}
    \item{n.alt}{The number of alternate alleles included in the test.}
    \item{n.sample.alt}{The number of samples with an observed alternate allele at any variant in the aggregate set.}
    
    If \code{burden.test} is \code{"BinomiRare"}:
    \item{n.carrier}{Number of individuals with at least one copy of the alternate allele}
    \item{n.D.carrier}{Number of cases with at least one copy of the alternate allele}
    \item{pval}{BinomiRare p-value}
    If \code{burden.test} is \code{"CMP"}:
    \item{n.carrier}{Number of individuals with at least one copy of the alternate allele}
    \item{n.D.carrier}{Number of cases with at least one copy of the alternate allele}
    \item{pval}{CMP p-value}
    \item{mid.pval}{CMP mid-p-value}

    \item{variantInfo}{A list with as many elements as aggregate tests performed. Each element of the list is a data.frame providing information on the variants used in the aggregate test with results presented in the corresponding row of \code{results}. Each of these data.frames has the following information:}
    \item{variant.id}{The variant ID}
    \item{chr}{The chromosome value}
    \item{pos}{The base pair position}
    \item{n.obs}{The number of samples with non-missing genotypes}
    \item{freq}{The estimated alternate allele frequency}
    \item{weight}{The weight assigned to the variant in the analysis.}
  


}

\author{Jiwon Lee, Tamar Sofer, Matthew P. Conomos, Stephanie M. Gogarten, Ken Rice, Chaoyu Yu, Han Chen}
\seealso{
     \code{\link{gen_id_list}} for conveniently creating an id list separated by group needed for \code{assocTestAggregate_split}
     \code{\link{fitNullModel}} for fitting the null mixed model needed as input to \code{assocTestSingle_split}.
     \code{\link{SeqVarIterator}} for creating the input object with genotypes.
     \code{\link{save_split_results} for saving results list as separate files by group}
}

\references{
Leal, S.M. & Li, B. (2008). Methods for Detecting Associations with Rare Variants for Common Diseases: Application to Analysis of Sequence Data. American Journal of Human Genetics, 83(3): 311-321.

Browning, S.R. & Madsen, B.E. (2009). A Groupwise Association Test for Rare Mutations Using a Weighted Sum Statistic. PLoS Genetics, 5(2): e1000384.

Wu, M.C, Lee, S., Cai, T., Li, Y., Boehnke, M., & Lin, X. (2011). Rare-Variant Association Testing for Sequencing Data with the Sequence Kernel Association Test. American Journal of Human Genetics, 89(1): 82-93.

Lee, S. et al. (2012). Optimal Unified Approach for Rare-Variant Association Testing with Application to Small-Sample Case-Control Whole-Exome Sequencing Studies. American Journal of Human Genetics, 91(2): 224-237.
}

\examples{
library(SeqVarTools)
library(Biobase)
library(GenomicRanges)

# open a sequencing GDS file
gdsfile <- seqExampleFileName("gds")
gds <- seqOpen(gdsfile)

# simulate some phenotype data
data(pedigree)
pedigree <- pedigree[match(seqGetData(gds, "sample.id"), pedigree$sample.id),]
pedigree$outcome <- rbinom(nrow(pedigree), 1, 0.4)

# Create id_list
id_list <- gen_id_list(pedigree, id_var="sample.id", group_var="sex")


# construct a SeqVarData object
seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(pedigree))

# fit the null model
nullmod <- fitNullModel(seqData, outcome="outcome", covars="sex", family="binomial")

# Range Iterator
gr <- GRanges(seqnames=rep(1,3), ranges=IRanges(start=c(1e6, 2e6, 3e6), width=1e6))
iterator <- SeqVarRangeIterator(seqData, variantRanges=gr)
assoc <- assocTestAggregate_split(iterator, nullmod, id_list=id_list, burden.test="BinomiRare", AF.max=0.01)

# Window Iterator
seqSetFilterChrom(seqData, include="22")
iterator <- SeqVarWindowIterator(seqData)
assoc <- assocTestAggregate_split(iterator, nullmod, id_list=id_list, burden.test="BinomiRare", AF.max=0.01)

# save results for each group separately
save_split_results(assoc, output_prefix='test')


seqClose(seqData)
}
\keyword{association}
