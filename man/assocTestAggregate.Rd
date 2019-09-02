\name{assocTestAggregate}
\alias{assocTestAggregate}
\alias{assocTestAggregate-methods}
\alias{assocTestAggregate,SeqVarIterator-method}
\alias{assocTestAggregate,GenotypeIterator-method}

\title{Aggregate Association Testing}
\description{\code{assocTestAggregate} performs aggregate association tests using the null model fit with \code{\link{fitNullModel}}.}

\usage{
\S4method{assocTestAggregate}{SeqVarIterator}(gdsobj, null.model, AF.max=1,
                   weight.beta=c(1,1), weight.user=NULL,
                   test=c("Burden", "SKAT", "fastSKAT", "SMMAT", "SKATO"),
                   burden.test=c("Score", "Wald"),
                   neig = 200, ntrace = 500,
                   rho = seq(from = 0, to = 1, by = 0.1),
                   sparse=TRUE, imputed=FALSE,
                   male.diploid=TRUE, genome.build=c("hg19", "hg38"),
                   verbose=TRUE)
\S4method{assocTestAggregate}{GenotypeIterator}(gdsobj, null.model, AF.max=1,
                   weight.beta=c(1,1), weight.user=NULL,
                   test=c("Burden", "SKAT", "fastSKAT", "SMMAT", "SKATO"),
                   burden.test=c("Score", "Wald"),
                   neig = 200, ntrace = 500,
                   rho = seq(from = 0, to = 1, by = 0.1),
                   male.diploid=TRUE, verbose=TRUE)
}

\arguments{
    \item{gdsobj}{An object of class \code{\link{SeqVarIterator}} from the package \pkg{\link{SeqVarTools}} containing the genotype data for the variants and samples to be used for the analysis.}
    \item{null.model}{A null model object returned by \code{fitNullModel}.}
    \item{AF.max}{A numeric value specifying the upper bound on the alternate allele frequency for variants to be included in the analysis.}
    \item{weight.beta}{A numeric vector of length two specifying the two parameters of the Beta distribution used to determine variant weights; weights are given by \code{dbeta(MAF, a, b)}, where MAF is the minor allele frequency, and a and b are the two parameters specified here. \code{weight.beta = c(1,25)} gives the Wu weights; \code{weight.beta = c(0.5, 0.5)} is proportional to the Madsen-Browning weights; and \code{weight.beta = c(1,1)} gives a weight of 1 to all variants. This input is ignored when \code{weight.user} is not \code{NULL}.}
    \item{weight.user}{A character string specifying the name of a variable to be used as variant weights. This variable can be in either 1) the variantData slot of \code{gdsobj} or 2) the \code{\link{mcols}} of the \code{\link{GRanges}} or \code{\link{GRangesList}} object used to create \code{gdsobj} (when \code{gdsobj} is a \code{link{SeqVarRangeIterator}} or \code{link{SeqVarListIterator}}).  When left \code{NULL} (the default), the weights specified by \code{weight.beta} will be used.}
    \item{test}{A character string specifying the type of test to be performed. The possibilities are \code{"Burden"} (default), \code{"SKAT"}, \code{"fastSKAT"}, \code{"SMMAT"}, or \code{"SKATO"}.}
    \item{burden.test}{A character string specifying the type of Burden test to perform when \code{test = "Burden"}.  The possibilities are \code{"Score"},  \code{"Wald"}, \code{"BinomiRare"}, and \code{"CMP"}. \code{"Score"} can be used for any \code{null.model}. \code{"Wald"} can not be used when the \code{null.model} is from a mixed model with a binary outcome variable.}
    %\item{burden.test}{A character string specifying the type of Burden test to perform when \code{test} = "Burden".  The possibilities are "Score", "Wald", and "Firth". "Score" can be used for any \code{null.model}. "Wald" can not be used when the \code{null.model} is from a mixed model with a binary outcome variable. "Firth" can only be used when the \code{null.model} is from a logistic regression with a binary outcome variable.}
    \item{neig}{The number eigenvalues to approximate by using random projections for calculating p-values with fastSKAT; default is 200. See 'Details' for more information.}
    \item{ntrace}{The number of vectors to sample when using random projections to estimate the trace needed for p-value calculation with fastSKAT; default is 500. See 'Details' for more information.}
    \item{rho}{A numeric value (or vector of numeric values) in \code{[0,1]} specifying the rho parameter when using \code{test == "SKATO"}; these are the values for which SKAT-O is performed, defining the search space for the optimal \code{rho}. If \code{rho = 0}, this is equivalent to a standard SKAT test; if \code{rho = 1}, this is equivalent to a score burden test.}
    %\item{pval.method}{A character string specifying which method to use to calculate SKAT p-values. \code{"davies"} (the default) uses numerical integration; \code{"kuonen"} uses a saddlepoint method; and \code{"liu"} uses a moment matching approximation. If the davies method generates an error, kuonen is tried, and then liu as a last resort.}
    \item{sparse}{Logical indicator of whether to read genotypes as sparse Matrix objects; the default is \code{TRUE}. Set this to \code{FALSE} if the alternate allele dosage of the genotypes in the test are not expected to be mostly 0.}
    \item{imputed}{Logical indicator of whether to read dosages from the "DS" field containing imputed dosages instead of counting the number of alternate alleles.}
  \item{male.diploid}{Logical for whether males on sex chromosomes are coded as diploid.}
  \item{genome.build}{A character sting indicating genome build; used to identify pseudoautosomal regions on the X and Y chromosomes.}
    \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is \code{TRUE}.}
}

\details{
    The type of aggregate unit tested depends on the class of iterator
    used for \code{gdsobj}. Options include sliding windows, specific
    ranges of variants or selection of individual variants (ranges with
    width 1). See \code{\link{SeqVarIterator}} for more details.
  
    The effect size estimate is for each copy of the alternate allele.
    For multiallelic variants, each alternate allele is tested separately.

    Monomorphic variants (including variants where every sample is a heterozygote) are always omitted from the aggregate unit prior to testing.
  
    Somewhat similarly to SKAT-O, the variant Set Mixed Model Association Test (SMMAT, Chen et al., 2019) combines the burden test p-value with an adjusted SKAT (which is asymptotically independent of the burden test) p-value using a chi-square distribution with 4df from Fisher's method.

    SKAT and SMMAT will attempt to use Davies' method (i.e. integration) to calculate p-values; if an error occurs in integration or the reported p-values are too small that they are unreliable (i.e. near machine epsilon), then the saddlepoint approximation will instead be used to calculate the p-values.

    The fastSKAT method of Lumley et al. (2018) uses random matrix theory to speed up the computation of SKAT p-values. When \code{test = "fastSKAT"}, the function attempts to inteligently determine which p-value calculation approach to use for each aggregation unit: (1) if min(number samples, number variants) is small enough, then the standard SKAT p-value calculation is used; (2) if min(number samples, number variants) is too large for standard SKAT, but small enough to explicitly compute the genotype covariance matrix, random projections are used to approximate the eigenvalues of the covariance matrix, and the fastSKAT p-value calculation is used; (3) if min(number samples, number variants) is too big to explicitly compute the genotype covariance matrix, random projections are used to approximate both the eigenvalues and the trace of the covariance matrix, and the fastSKAT p-value calculation is used.)
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
    If \code{test} is \code{"Burden"}:
    %% \item{burden.skew}{The skewness of the burden value for all samples.}
    
    If \code{burden.test} is "Score":
    \item{Score}{The value of the score function}
    \item{Score.SE}{The estimated standard error of the Score}
    \item{Score.Stat}{The score Z test statistic}
    \item{Score.pval}{The score p-value}
    If \code{burden.test} is \code{"Wald"}:
    \item{Est}{The effect size estimate for a one unit increase in the burden value}
    \item{Est.SE}{The estimated standard error of the effect size estimate}    
    \item{Wald.Stat}{The Wald Z test statistic}
    \item{Wald.pval}{The Wald p-value}
     If \code{burden.test} is \code{"BinomiRare"}:
    \item{n.carrier}{Number of individuals with at least one copy of the alternate allele}
    \item{n.D.carrier}{Number of cases with at least one copy of the alternate allele}
    \item{pval}{BinomiRare p-value}
    If \code{burden.test} is \code{"CMP"}:
    \item{n.carrier}{Number of individuals with at least one copy of the alternate allele}
    \item{n.D.carrier}{Number of cases with at least one copy of the alternate allele}
    \item{pval}{CMP p-value}
    \item{mid.pval}{CMP mid-p-value}
    %% If \code{burden.test} is "Firth":
    %% \item{Est}{The effect size estimate for a one unit increase in the burden value}
    %% \item{SE}{The estimated standard error of the effect size estimate}
    %% \item{Firth.stat}{The Firth test statistic}
    %% \item{Firth.pval}{The Firth p-value}
    If \code{test} is \code{"SKAT"} or \code{"fastSKAT"}:
    \item{Q}{The SKAT test statistic.}
    \item{pval}{The SKAT p-value.}
    \item{err}{Takes value 1 if there was an error in calculating the p-value; takes the value 2 when multiple random projections were required to get a good approximation from fastSKAT (the reported p-value is likely still reliable); 0 otherwise.}
    \item{pval.method}{The p-value calculation method used. When standard SKAT is used, one of "integration" or "saddlepoint"; when fastSKAT random projections are used to approximate eigenvalues of the genotype covariance matrix, one of "ssvd_integration" or "ssvd_saddlepoint"; when fastSKAT random projections are used to approximate both the eigenvalues and the trace of the genotype covariance matrix, one of "rsvd_integration" or "rsvd_saddlepoint".}

    If \code{test} is \code{"SMMAT"}:
    \item{pval_burden}{The burden test p-value}
    \item{pval_theta}{The p-value of the adjusted SKAT test (which is asymptotically independent of the burden test)}
    \item{pval_SMMAT}{The SMMAT p-value after combining pval_burden and pval_theta using Fisher's method.}
    \item{err}{Takes value 1 if there was an error calculating the SMMAT p-value; 0 otherwise. If \code{err=1}, \code{pval_SMMAT} is set to \code{pval_burden}.}

    If \code{test} is \code{"SKATO"}:
    \item{Q_rho}{The SKAT test statistic for the value of rho specified. There will be as many of these variables as there are rho values chosen.}
    \item{pval_rho}{The SKAT p-value for the value of rho specified.  There will be as many of these variables as there are rho values chosen.}
    \item{err_rho}{Takes value 1 if there was an error in calculating the p-value for the value of rho specified when using the "kuonen" or "davies" methods; 0 otherwise. When there is an error, the p-value returned is from the "liu" method. There will be as many of these variables as there are rho values chosen.}
    \item{min.pval}{The minimum p-value among the p-values calculated for each choice of rho.}
    \item{opt.rho}{The optimal rho value; i.e. the rho value that gave the minimum p-value.}
    \item{pval_SKATO}{The SKAT-O p-value after adjustment for searching across multiple rho values.}

    \item{variantInfo}{A list with as many elements as aggregate tests performed. Each element of the list is a data.frame providing information on the variants used in the aggregate test with results presented in the corresponding row of \code{results}. Each of these data.frames has the following information:}
    \item{variant.id}{The variant ID}
    \item{chr}{The chromosome value}
    \item{pos}{The base pair position}
    \item{allele.index}{The index of the alternate allele. For biallelic variants, this will always be 1.}
    \item{n.obs}{The number of samples with non-missing genotypes}
    \item{freq}{The estimated alternate allele frequency}
    \item{MAC}{The minor allele count. For multiallelic variants, "minor" is determined by comparing the count of the alternate allele specified by \code{allele.index} with the sum of all other alleles.}
    \item{weight}{The weight assigned to the variant in the analysis.}
}

\author{Matthew P. Conomos, Stephanie M. Gogarten, Thomas Lumley, Tamar Sofer, Ken Rice, Chaoyu Yu, Han Chen}

\references{
Leal, S.M. & Li, B. (2008). Methods for Detecting Associations with Rare Variants for Common Diseases: Application to Analysis of Sequence Data. American Journal of Human Genetics, 83(3): 311-321.

Browning, S.R. & Madsen, B.E. (2009). A Groupwise Association Test for Rare Mutations Using a Weighted Sum Statistic. PLoS Genetics, 5(2): e1000384.

Wu, M.C, Lee, S., Cai, T., Li, Y., Boehnke, M., & Lin, X. (2011). Rare-Variant Association Testing for Sequencing Data with the Sequence Kernel Association Test. American Journal of Human Genetics, 89(1): 82-93.

Lee, S. et al. (2012). Optimal Unified Approach for Rare-Variant Association Testing with Application to Small-Sample Case-Control Whole-Exome Sequencing Studies. American Journal of Human Genetics, 91(2): 224-237.

Chen, H., Huffman, J. E., Brody, J. A., Wang, C., Lee, S., Li, Z., ... & Blangero, J. (2019). Efficient variant set mixed model association tests for continuous and binary traits in large-scale whole-genome sequencing studies. The American Journal of Human Genetics, 104(2), 260-274.

Lumley, T., Brody, J., Peloso, G., Morrison, A., & Rice, K. (2018). FastSKAT: Sequence kernel association tests for very large sets of markers. Genetic epidemiology, 42(6), 516-527.
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
pedigree$outcome <- rnorm(nrow(pedigree))

# construct a SeqVarData object
seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(pedigree))

# fit the null model
nullmod <- fitNullModel(seqData, outcome="outcome", covars="sex")

# burden test - Range Iterator
gr <- GRanges(seqnames=rep(1,3), ranges=IRanges(start=c(1e6, 2e6, 3e6), width=1e6))
iterator <- SeqVarRangeIterator(seqData, variantRanges=gr)
assoc <- assocTestAggregate(iterator, nullmod, test="Burden")
assoc$results
lapply(assoc$variantInfo, head)

# SKAT test - Window Iterator
seqSetFilterChrom(seqData, include="22")
iterator <- SeqVarWindowIterator(seqData)
assoc <- assocTestAggregate(iterator, nullmod, test="SKAT")
head(assoc$results)
head(assoc$variantInfo)

# SKAT-O test - List Iterator
seqResetFilter(iterator)
gr <- GRangesList(
  GRanges(seqnames=rep(22,2), ranges=IRanges(start=c(16e6, 17e6), width=1e6)),
  GRanges(seqnames=rep(22,2), ranges=IRanges(start=c(18e6, 20e6), width=1e6)))
iterator <- SeqVarListIterator(seqData, variantRanges=gr)
assoc <- assocTestAggregate(iterator, nullmod, test="SKAT", rho=seq(0, 1, 0.25))
assoc$results
assoc$variantInfo

# user-specified weights - option 1
seqResetFilter(iterator)
variant.id <- seqGetData(gds, "variant.id")
weights <- data.frame(variant.id, weight=runif(length(variant.id)))
variantData(seqData) <- AnnotatedDataFrame(weights)
iterator <- SeqVarListIterator(seqData, variantRanges=gr)
assoc <- assocTestAggregate(iterator, nullmod, test="Burden", weight.user="weight")
assoc$results
assoc$variantInfo

# user-specified weights - option 2
seqResetFilter(iterator)
variantData(seqData)$weight <- NULL
gr <- GRangesList(
  GRanges(seqnames=rep(22,2), ranges=IRanges(start=c(16e6, 17e6), width=1e6), weight=runif(2)),
  GRanges(seqnames=rep(22,2), ranges=IRanges(start=c(18e6, 20e6), width=1e6), weight=runif(2)))
iterator <- SeqVarListIterator(seqData, variantRanges=gr)
assoc <- assocTestAggregate(iterator, nullmod, test="Burden", weight.user="weight")
assoc$results
assoc$variantInfo

seqClose(seqData)
}
\keyword{association}
