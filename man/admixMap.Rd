\name{admixMap}
\alias{admixMap}

\title{admixMap}
\description{
Run admixture analyses
}
\usage{
admixMap(admixDataList, null.model, male.diploid=TRUE,
         genome.build=c("hg19", "hg38"),
         BPPARAM=bpparam(), verbose=TRUE)
}

\arguments{
  \item{admixDataList}{named list of \code{\link{GenotypeIterator}} or \code{\link{SeqVarIterator}} objects for each ancestry}
    \item{null.model}{A null model object returned by \code{\link{fitNullModel}}.}
    \item{male.diploid}{Logical for whether males on sex chromosomes are coded as diploid. Default is `male.diploid=TRUE`, meaning sex chromosome genotypes for males have values 0/2. If the input object codes males as 0/1 on sex chromosomes, set `male.diploid=FALSE`.}
  \item{genome.build}{A character sting indicating genome build; used to identify pseudoautosomal regions on the X and Y chromosomes. These regions are not treated as sex chromosomes when calculating allele frequencies.}
  \item{BPPARAM}{A \code{\link{BiocParallelParam}} object to process blocks of variants in parallel. If not provided, the default back-end returned by \code{\link{bpparam}} will be used.}
    \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is TRUE.}
}

\details{
This function is used with local ancestry results such as those obtained from RFMix. RFMix output may be converted to PLINK format, and then to GDS with \code{\link[SNPRelate]{snpgdsBED2GDS}}.

\code{admixDataList} should have one value for each ancestry to be included in the test. The sum of all local ancestries at a particular locus must add up to 2, so if there are K ancestry groups, then only K-1 genotypes can be included since one local ancestry count can be written as a linear combination of all of the other local ancestry counts, resulting in collinearity and a matrix that won't be invertible.

See the example for how one might set up the \code{admixDataList} object. List names will propagate to the output file.

    \code{admixMap} uses the \code{\link{BiocParallel}} package to process iterator chunks in parallel. See the \code{\link{BiocParallel}} documentation for more information on the default behaviour of \code{\link{bpparam}} and how to register different parallel backends. If serial execution is desired, set \code{BPPARAM=BiocParallel::SerialParam()}. Note that parallel execution requires more RAM than serial execution.
}

\value{A data.frame where each row refers to a different variant with the columns:
    \item{variant.id}{The variant ID}
    \item{chr}{The chromosome value}
    \item{pos}{The base pair position}
    \item{n.obs}{The number of samples with non-missing genotypes}
    \item{*.freq}{The estimated frequency of alleles derived from each ancestry at that variant}
    \item{*.Est}{The effect size estimate for each additional copy of an allele derived from each ancestry, relative to the reference ancestry}
    \item{*.SE}{The estimated standard error of the effect size estimate for each ancestry}
    \item{Joint.Stat}{The chi-square Wald test statistic for the joint test of all local ancestry terms}
    \item{Joint.pval}{The Wald p-value for the joint test of all local ancestry terms}
}

p-values that are calculated using \code{pchisq} and are smaller than \code{.Machine\$double.xmin} are set to \code{.Machine\$double.xmin}.

\author{Matthew P. Conomos, Lisa Brown, Stephanie M. Gogarten, Tamar Sofer, Ken Rice, Chaoyu Yu}

\seealso{\code{\link{GenotypeIterator}}, \code{\link{fitNullModel}}, \code{\link{assocTestSingle}}}

\examples{
library(GWASTools)

# option 1: one GDS file per ancestry
afrfile <- system.file("extdata", "HapMap_ASW_MXL_local_afr.gds", package="GENESIS")
amerfile <- system.file("extdata", "HapMap_ASW_MXL_local_amer.gds", package="GENESIS")
eurfile <- system.file("extdata", "HapMap_ASW_MXL_local_eur.gds", package="GENESIS")
files <- list(afr=afrfile, amer=amerfile, eur=eurfile)
gdsList <- lapply(files, GdsGenotypeReader)

# make ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(data.frame(
    scanID=getScanID(gdsList[[1]]), stringsAsFactors=FALSE))

# generate a phenotype
set.seed(4)
nsamp <- nrow(scanAnnot)
scanAnnot$pheno <- rnorm(nsamp, mean=0, sd=1)
set.seed(5)
scanAnnot$covar <- sample(0:1, nsamp, replace=TRUE)

genoDataList <- lapply(gdsList, GenotypeData, scanAnnot=scanAnnot)

# iterators
# if we have 3 ancestries total, only 2 should be included in test
genoIterators <- lapply(genoDataList[1:2], GenotypeBlockIterator)

# fit the null mixed model
null.model <- fitNullModel(scanAnnot, outcome="pheno", covars="covar")

# run the association test
myassoc <- admixMap(genoIterators, null.model,
                    BPPARAM=BiocParallel::SerialParam())
head(myassoc)

lapply(genoDataList, close)


# option 2: create a single file with multiple ancestries
# first, get dosages from all ancestries
library(gdsfmt)
dosages <- lapply(files, function(f) {
    gds <- openfn.gds(f)
    geno <- read.gdsn(index.gdsn(gds, "genotype"))
    closefn.gds(gds)
    geno
})
lapply(dosages, dim)

# create a new file with three dosage matrices, keeping all
# sample and snp nodes from one original file
tmpfile <- tempfile()
file.copy(afrfile, tmpfile)
gds <- openfn.gds(tmpfile, readonly=FALSE)
delete.gdsn(index.gdsn(gds, "genotype"))
add.gdsn(gds, "dosage_afr", dosages[["afr"]])
add.gdsn(gds, "dosage_amer", dosages[["amer"]])
add.gdsn(gds, "dosage_eur", dosages[["eur"]])
closefn.gds(gds)
cleanup.gds(tmpfile)

# read in GDS data, specifying the node for each ancestry
gds <- openfn.gds(tmpfile)
gds
genoDataList <- list()
for (anc in c("afr", "amer", "eur")){
  gdsr <- GdsGenotypeReader(gds, genotypeVar=paste0("dosage_", anc))
  genoDataList[[anc]] <- GenotypeData(gdsr, scanAnnot=scanAnnot)
}

# iterators
genoIterators <- lapply(genoDataList[1:2], GenotypeBlockIterator)

# run the association test
myassoc <- admixMap(genoIterators, null.model,
                    BPPARAM=BiocParallel::SerialParam())

close(genoDataList[[1]])
unlink(tmpfile)
}

\references{
Brown, L.A. et al. (2017). Admixture Mapping Identifies an Amerindian Ancestry Locus Associated with Albuminuria in Hispanics in the United States. J Am Soc Nephrol. 28(7):2211-2220.

Maples, B.K. et al. (2013). RFMix: a discriminative modeling approach for rapid and robust local-ancestry inference. Am J Hum Genet. 93(2):278-88.
}

\keyword{mixed model}
\keyword{association}
