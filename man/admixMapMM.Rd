\name{admixMapMM}
\alias{admixMapMM}

\title{admixMapMM}
\description{
Run admixture analyses
}
\usage{
admixMapMM(admixDataList, nullMMobj, snp.include = NULL,
           chromosome = NULL, snp.block.size = 5000,
           verbose = TRUE)
}

\arguments{
  \item{admixDataList}{named list of \code{\link{GenotypeData}} objects for each ancestry}
    \item{nullMMobj}{A null model object returned by \code{fitNullMM}.}
    \item{snp.include}{A vector of SNP IDs to include in the analysis.  If NULL, see \code{chromosome} for further details.}
    \item{chromosome}{A vector of integers specifying which chromosomes to analyze.  This parameter is only considered when \code{snp.include} is NULL; if \code{chromosome} is also NULL, then all SNPs are included.}
    \item{snp.block.size}{The number of SNPs to read-in/analyze at once.  The default value is 5000.}
    \item{verbose}{Logical indicator of whether updates from the function should be printed to the console; the default is TRUE.}
}

\details{
\code{admixDataList} should have one value for each ancestry. See the example for how one might set up this object. List names will propagate to the output file.
}

\value{data frame with admixture mapping results}

\author{Matt Conomos}

\seealso{\code{\link{GenotypeData}}, \code{\link{fitNullMM}}, \code{\link{assocTestMM}}}

\examples{
library(GWASTools)
library(gdsfmt)

# create file with multiple ancestries
gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
tmpfile <- tempfile()
file.copy(gdsfile, tmpfile)
gds <- openfn.gds(tmpfile, readonly=FALSE)
nsnp <- objdesp.gdsn(index.gdsn(gds, "snp.id"))$dim
nsamp <- objdesp.gdsn(index.gdsn(gds, "sample.id"))$dim
dosage_eur <- sample(0:2, nsnp*nsamp, replace=TRUE)
dosage_afr <- ifelse(dosage_eur == 2, 0, sample(0:1, nsnp*nsamp, replace=TRUE))
dosage_amer <- 2 - dosage_eur - dosage_afr
add.gdsn(gds, "dosage_eur", matrix(dosage_eur, nrow=nsamp, ncol=nsnp))
add.gdsn(gds, "dosage_afr", matrix(dosage_afr, nrow=nsamp, ncol=nsnp))
add.gdsn(gds, "dosage_amer", matrix(dosage_amer, nrow=nsamp, ncol=nsnp))
closefn.gds(gds)
        
# read GRM
pcrfile <- system.file("extdata", "HapMap_ASW_MXL_pcrelate.gds", package="GENESIS")
pcr <- openfn.gds(pcrfile)
mypcrel <- pcrelateMakeGRM(pcr)
closefn.gds(pcr)

# generate a phenotype
set.seed(4)
pheno <- rnorm(nsamp, mean = 0, sd = 1)
covar <- sample(0:1, nsamp, replace=TRUE)

# make ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID = rownames(mypcrel), 
              covar, pheno, stringsAsFactors=FALSE))

# read in GDS data
gds <- openfn.gds(tmpfile)
genoDataList <- list()
for (anc in c("eur", "afr", "amer")){
  gdsr <- GdsGenotypeReader(gds, genotypeVar=paste0("dosage_", anc))
  genoDataList[[anc]] <- GenotypeData(gdsr, scanAnnot=scanAnnot)
}
     
# fit the null mixed model
nullmod <- fitNullMM(scanData = scanAnnot, outcome = "pheno", covars = "covar", covMatList = mypcrel)

# run the association test
myassoc <- admixMapMM(genoDataList, nullMMobj = nullmod)

close(genoDataList[[1]])
unlink(tmpfile)
}

\keyword{mixed model}
\keyword{association}
