## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE)

## ---- echo=FALSE, results='hide'----------------------------------------------
library(GENESIS)
library(GWASTools)

# file path to GDS file
gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
# read in GDS data
HapMap_geno <- GdsGenotypeReader(filename = gdsfile)
# create a GenotypeData class object
HapMap_genoData <- GenotypeData(HapMap_geno)
# load saved matrix of KING-robust estimates
data("HapMap_ASW_MXL_KINGmat")

# run PC-AiR
mypcair <- pcair(HapMap_genoData,
                 kinobj = HapMap_ASW_MXL_KINGmat,
                 divobj = HapMap_ASW_MXL_KINGmat,
                 verbose = FALSE)
mypcs <- mypcair$vectors[,1,drop=FALSE]

# create a GenotypeBlockIterator object
HapMap_genoData <- GenotypeBlockIterator(HapMap_genoData)
# run PC-Relate
mypcrel <- pcrelate(HapMap_genoData, pcs = mypcs,
                    training.set = mypcair$unrels,
                    verbose = FALSE)

# generate a phenotype
set.seed(4)
pheno <- 0.2*mypcs + rnorm(mypcair$nsamp, mean = 0, sd = 1)

## -----------------------------------------------------------------------------
# mypcair contains PCs from a previous PC-AiR analysis
# pheno is a vector of Phenotype values

# make a data.frame
mydat <- data.frame(scanID = mypcair$sample.id, pc1 = mypcair$vectors[,1],
                    pheno = pheno)
head(mydat)

# make ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(mydat)
scanAnnot

## ---- eval=FALSE--------------------------------------------------------------
#  geno <- MatrixGenotypeReader(genotype = genotype, snpID = snpID,
#                               chromosome = chromosome, position = position,
#                               scanID = scanID)
#  genoData <- GenotypeData(geno)

## ---- eval=FALSE--------------------------------------------------------------
#  geno <- GdsGenotypeReader(filename = "genotype.gds")
#  genoData <- GenotypeData(geno)

## ---- eval=FALSE--------------------------------------------------------------
#  snpgdsBED2GDS(bed.fn = "genotype.bed",
#                bim.fn = "genotype.bim",
#                fam.fn = "genotype.fam",
#                out.gdsfn = "genotype.gds")

## ---- eval=FALSE--------------------------------------------------------------
#  # read in GDS data
#  gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")
#  #HapMap_geno <- GdsGenotypeReader(filename = gdsfile)

## -----------------------------------------------------------------------------
# create a GenotypeData class object with paired ScanAnnotationDataFrame
HapMap_genoData <- GenotypeData(HapMap_geno, scanAnnot = scanAnnot)
HapMap_genoData

## -----------------------------------------------------------------------------
# mypcrel contains Kinship Estimates from a previous PC-Relate analysis
myGRM <- pcrelateToMatrix(mypcrel)
myGRM[1:5,1:5]

## -----------------------------------------------------------------------------
# fit the null mixed model
nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = "pc1",
                        cov.mat = myGRM, family = gaussian)

## ---- eval=FALSE--------------------------------------------------------------
#  nullmod <- fitNullModel(scanAnnot, outcome = "pheno",
#                          covars = c("pc1","pc2","sex","age"),
#                          cov.mat = myGRM, family = gaussian)

## ---- eval=FALSE--------------------------------------------------------------
#  nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = "pc1",
#                          cov.mat = list("GRM" = myGRM, "House" = H),
#                          family = gaussian)

## ---- eval=FALSE--------------------------------------------------------------
#  nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = "pc1",
#                          cov.mat = myGRM, family = gaussian,
#                          group.var = "study")

## ---- eval=FALSE--------------------------------------------------------------
#  nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = "pc1",
#                          cov.mat = myGRM, family = binomial)

## -----------------------------------------------------------------------------
genoIterator <- GenotypeBlockIterator(HapMap_genoData, snpBlock=5000)

## -----------------------------------------------------------------------------
assoc <- assocTestSingle(genoIterator, null.model = nullmod)

## ---- eval = FALSE------------------------------------------------------------
#  # mysnps is a vector of snpID values for the SNPs we want to test
#  genoIterator <- GenotypeBlockIterator(HapMap_genoData, snpInclude=mysnps)
#  assoc <- assocTestSingle(genoIterator, null.model = nullmod)

## -----------------------------------------------------------------------------
head(assoc)

