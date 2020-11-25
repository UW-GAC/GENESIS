## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, fig.height = 6, fig.width = 6)

## ----vcf2gds------------------------------------------------------------------
library(SeqArray)
vcffile <- system.file("extdata", "1KG", 
                       paste0("1KG_phase3_subset_chr", 1:22, ".vcf.gz"), 
                       package="GENESIS")
gdsfile <- tempfile()
seqVCF2GDS(vcffile, gdsfile, verbose=FALSE)
gds <- seqOpen(gdsfile)
gds

## ----seqvardata---------------------------------------------------------------
library(GENESIS)
library(Biobase)
library(SeqVarTools)

data(sample_annotation_1KG)
annot <- sample_annotation_1KG
head(annot)

# simulate some phenotype data
set.seed(4)
annot$outcome <- rnorm(nrow(annot))
metadata <- data.frame(labelDescription=c("sample id", 
                                          "1000 genomes population", 
                                          "sex", 
                                          "simulated phenotype"),
                       row.names=names(annot))
annot <- AnnotatedDataFrame(annot, metadata)

all.equal(annot$sample.id, seqGetData(gds, "sample.id"))
seqData <- SeqVarData(gds, sampleData=annot)

## ----seed, include=FALSE------------------------------------------------------
# set seed for LD pruning
set.seed(100)

## ----king---------------------------------------------------------------------
library(SNPRelate)

# LD pruning to get variant set
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, 
                          ld.threshold=sqrt(0.1), verbose=FALSE)
pruned <- unlist(snpset, use.names=FALSE)

king <- snpgdsIBDKING(gds, snp.id=pruned, verbose=FALSE)
kingMat <- king$kinship
dimnames(kingMat) <- list(king$sample.id, king$sample.id)

## ----pcair--------------------------------------------------------------------
pcs <- pcair(seqData, 
             kinobj=kingMat, kin.thresh=2^(-9/2),
             divobj=kingMat, div.thresh=-2^(-9/2),
             snp.include=pruned)

## ----pcair_plot, fig.width=7, out.width="100%"--------------------------------
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(GGally)

pc.df <- as.data.frame(pcs$vectors)
names(pc.df) <- paste0("PC", 1:ncol(pcs$vectors))
pc.df$sample.id <- row.names(pcs$vectors)
pc.df <- left_join(pc.df, pData(annot), by="sample.id")

pop.cols <- setNames(brewer.pal(12, "Paired"),
    c("ACB", "ASW", "CEU", "GBR", "CHB", "JPT", 
      "CLM", "MXL", "LWK", "YRI", "GIH", "PUR"))

ggplot(pc.df, aes(PC1, PC2, color=Population)) + geom_point() +
    scale_color_manual(values=pop.cols)

## ----parcoord, fig.wide=TRUE, fig.height=4, fig.width=10----------------------
ggparcoord(pc.df, columns=1:10, groupColumn="Population", scale="uniminmax") +
    scale_color_manual(values=pop.cols) +
    xlab("PC") + ylab("")

## ----pcrelate-----------------------------------------------------------------
seqSetFilter(seqData, variant.id=pruned)
iterator <- SeqVarBlockIterator(seqData, variantBlock=20000, verbose=FALSE)
pcrel <- pcrelate(iterator, pcs=pcs$vectors[,1:2], training.set=pcs$unrels)
seqResetFilter(seqData, verbose=FALSE)

## ----pcrelate_plot------------------------------------------------------------
kinship <- pcrel$kinBtwn

ggplot(kinship, aes(k0, kin)) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
    geom_point(alpha=0.5) +
    ylab("kinship estimate") +
    theme_bw()

## ----null_model---------------------------------------------------------------
# add PCs to sample annotation in SeqVarData object
annot <- AnnotatedDataFrame(pc.df)
sampleData(seqData) <- annot

# covariance matrix from pcrelate output
grm <- pcrelateToMatrix(pcrel, scaleKin=2)

# fit the null model
nullmod <- fitNullModel(seqData, outcome="outcome", 
                        covars=c("sex", "Population", paste0("PC", 1:2)),
                        cov.mat=grm, verbose=FALSE)

## ----assoc_single-------------------------------------------------------------
# select chromosome 1
seqSetFilterChrom(seqData, include=1)

iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
assoc <- assocTestSingle(iterator, nullmod, verbose=FALSE)
head(assoc)

## ----assoc_single_qq----------------------------------------------------------
qqPlot <- function(pval) {
    pval <- pval[!is.na(pval)]
    n <- length(pval)
    x <- 1:n
    dat <- data.frame(obs=sort(pval),
                      exp=x/n,
                      upper=qbeta(0.025, x, rev(x)),
                      lower=qbeta(0.975, x, rev(x)))
    
    ggplot(dat, aes(-log10(exp), -log10(obs))) +
        geom_line(aes(-log10(exp), -log10(upper)), color="gray") +
        geom_line(aes(-log10(exp), -log10(lower)), color="gray") +
        geom_point() +
        geom_abline(intercept=0, slope=1, color="red") +
        xlab(expression(paste(-log[10], "(expected P)"))) +
        ylab(expression(paste(-log[10], "(observed P)"))) +
        theme_bw()
}    

qqPlot(assoc$Score.pval)

## -----------------------------------------------------------------------------
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# return the variants on chromosome 1 as a GRanges object
seqSetFilterChrom(seqData, include=1)
gr <- granges(gds)

# find variants that overlap with each gene
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gr <- renameSeqlevels(gr, paste0("chr", seqlevels(gr)))
ts <- transcriptsByOverlaps(txdb, gr, columns="GENEID")
# simplistic example - define genes as overlapping transcripts
genes <- reduce(ts)
genes <- renameSeqlevels(genes, sub("chr", "", seqlevels(genes)))

# create an iterator where each successive unit is a different gene
iterator <- SeqVarRangeIterator(seqData, variantRanges=genes, verbose=FALSE)

# do a burden test on the rare variants in each gene
assoc <- assocTestAggregate(iterator, nullmod, AF.max=0.05, test="Burden",
                            verbose=FALSE)

## -----------------------------------------------------------------------------
head(assoc$results)
head(assoc$variantInfo)

## ----close, echo=FALSE--------------------------------------------------------
seqClose(gds)
unlink(gdsfile)

