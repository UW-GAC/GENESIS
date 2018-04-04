\name{assocTestAggregate}
\alias{assocTestAggregate}
\alias{assocTestAggregate-methods}
\alias{assocTestAggregate,SeqVarIterator-method}

\title{Aggregate Association Testing}
\description{\code{assocTestAggregate} performs aggregate association tests using the null model fit with \code{\link{fitNullModel}}.}

\usage{
\S4method{assocTestAggregate}{SeqVarIterator}(gdsobj, null.model, AF.max=1,
                   weight.beta=c(1,1), weight.user=NULL,
                   test=c("Burden", "SKAT", "SMMAT"),
                   burden.test=c("Score", "Wald"), rho=0,
                   pval.method=c("davies", "kuonen", "liu"),
                   verbose=TRUE)
}

\arguments{
    \item{gdsobj}{An object of class \code{\link{SeqVarIterator}} from the package \pkg{\link{SeqVarTools}} containing the genotype data for the variants and samples to be used for the analysis.}
    \item{null.model}{A null model object returned by \code{fitNullModel}.}
    \item{AF.max}{A numeric value specifying the upper bound on the alternate allele frequency for variants to be included in the analysis.}
    \item{weight.beta}{A numeric vector of length two specifying the two parameters of the Beta distribution used to determine variant weights; weights are given by \code{dbeta(MAF, a, b)}, where MAF is the minor allele frequency, and a and b are the two parameters specified here. \code{weight.beta = c(1,25)} gives the Wu weights; \code{weight.beta = c(0.5, 0.5)} is proportional to the Madsen-Browning weights; and \code{weight.beta = c(1,1)} gives a weight of 1 to all variants. This input is ignored when \code{weight.user} is not \code{NULL}.}
    \item{weight.user}{A character string specifying the name of a variable in the variantData slot of \code{gdsobj} to be used as variant weights.  When left \code{NULL} (the default), the weights specified by \code{weight.beta} will be used.}
    \item{test}{A character string specifying the type of test to be performed. The possibilities are \code{"Burden"} (default), \code{"SKAT"}, or \code{"SMMAT"}.  When this is set to "SKAT" and the parameter \code{rho} has multiple values, a SKAT-O test is performed.}
    \item{burden.test}{A character string specifying the type of Burden test to perform when \code{test = "Burden"}.  The possibilities are \code{"Score"} and \code{"Wald"}. \code{"Score"} can be used for any \code{null.model}. \code{"Wald"} can not be used when the \code{null.model} is from a mixed model with a binary outcome variable.}
    %\item{burden.test}{A character string specifying the type of Burden test to perform when \code{test} = "Burden".  The possibilities are "Score", "Wald", and "Firth". "Score" can be used for any \code{null.model}. "Wald" can not be used when the \code{null.model} is from a mixed model with a binary outcome variable. "Firth" can only be used when the \code{null.model} is from a logistic regression with a binary outcome variable.}
    \item{rho}{A numeric value (or vector of numeric values) in \code{[0,1]} specifying the rho parameter for SKAT. When \code{rho = 0}, a standard SKAT test is performed. When \code{rho = 1}, a score burden test is performed. When \code{rho} is a vector of values, SKAT-O is performed using each of those values as the search space for the optimal \code{rho}.}
    \item{pval.method}{A character string specifying which method to use to calculate SKAT p-values. \code{"davies"} (the default) uses numerical integration; \code{"kuonen"} uses a saddlepoint method; and \code{"liu"} uses a moment matching approximation. If the davies method generates an error, kuonen is tried, and then liu as a last resort.}
    \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is \code{TRUE}.}
}

\details{
    The type of aggregate unit tested depends on the class of iterator used for \code{gdsobj}. Options include sliding windows, specific ranges of variants or selection of individual variants (ranges with width 1). See \code{\link{SeqVarIterator}} for more details.
  
    The effect size estimate is for each copy of the alternate allele.
    For multiallelic variants, each alternate allele is tested separately.

    The SMMAT test is a hybrid of SKAT and the burdren test.
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
    %% If \code{burden.test} is "Firth":
    %% \item{Est}{The effect size estimate for a one unit increase in the burden value}
    %% \item{SE}{The estimated standard error of the effect size estimate}
    %% \item{Firth.stat}{The Firth test statistic}
    %% \item{Firth.pval}{The Firth p-value}
    If \code{test} is \code{"SKAT"}:
    \item{Q_rho}{The SKAT test statistic for the value of rho specified. There will be as many of these variables as there are rho values chosen.}
    \item{pval_rho}{The SKAT p-value for the value of rho specified.  There will be as many of these variables as there are rho values chosen.}
    \item{err_rho}{Takes value 1 if there was an error in calculating the p-value for the value of rho specified when using the "kuonen" or "davies" methods; 0 otherwise. When there is an error, the p-value returned is from the "liu" method. There will be as many of these variables as there are rho values chosen.}
    When \code{length(rho) > 1} and SKAT-O is performed:
    \item{min.pval}{The minimum p-value among the p-values calculated for each choice of rho.}
    \item{opt.rho}{The optimal rho value; i.e. the rho value that gave the minimum p-value.}
    \item{pval_SKATO}{The SKAT-O p-value after adjustment for searching across multiple rho values.}

    If \code{test} is \code{"SMMAT"}:
    \item{pval_burden}{The burden test p-value}
    \item{pval_hybrid}{The SMMAT p-value}
    \item{err}{Takes value 1 if there was an error calculating the hybrid p-value; 0 otherwise. If \code{err=1}, \code{pval_hybrid} is set to \code{pval_burden}.}
    

    \item{variantInfo}{A list with as many elements as aggregate tests performed. Each element of the list is a data.frame providing information on the variants used in the aggregate test with results presented in the corresponding row of \code{results}. Each of these data.frames has the following information:}
    \item{variant.id}{The variant ID}
    \item{chr}{The chromosome value}
    \item{pos}{The base pair position}
    \item{n.obs}{The number of samples with non-missing genotypes}
    \item{freq}{The estimated alternate allele frequency}
    \item{weight}{The weight assigned to the variant in the analysis.}
}

\author{Matthew P. Conomos, Stephanie M. Gogarten, Tamar Sofer, Ken Rice, Chaoyu Yu}

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

# user-specified weights
seqResetFilter(iterator)
variant.id <- seqGetData(gds, "variant.id")
weights <- data.frame(variant.id, weight=runif(length(variant.id)))
variantData(seqData) <- AnnotatedDataFrame(weights)
iterator <- SeqVarListIterator(seqData, variantRanges=gr)
assoc <- assocTestAggregate(iterator, nullmod, test="Burden", weight.user="weight")
assoc$results
assoc$variantInfo

seqClose(seqData)
}
\keyword{association}
