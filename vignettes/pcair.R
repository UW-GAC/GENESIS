## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE)

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

## -----------------------------------------------------------------------------
gdsfile <- system.file("extdata", "HapMap_ASW_MXL_geno.gds", package="GENESIS")

## ----seed, include=FALSE------------------------------------------------------
# set seed for LD pruning
set.seed(100)

