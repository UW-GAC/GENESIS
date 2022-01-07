# GENESIS
GENetic EStimation and Inference in Structured samples (GENESIS): Statistical methods for analyzing genetic data from samples with population structure and/or relatedness

The GENESIS package provides methodology for estimating,
        inferring, and accounting for population and pedigree structure
        in genetic analyses.  The current implementation provides
        functions to perform PC-AiR (Conomos et al., 2015, Gen Epi) and PC-Relate 
        (Conomos et al., 2016, AJHG). PC-AiR performs a Principal Components 
        Analysis on genome-wide SNP data for the detection of population 
        structure in a sample that may contain known or cryptic relatedness. 
        Unlike standard PCA, PC-AiR accounts for relatedness in the sample 
        to provide accurate ancestry inference that is not confounded by 
        family structure. PC-Relate uses ancestry representative principal 
        components to adjust for population structure/ancestry and accurately 
        estimate measures of recent genetic relatedness such as kinship 
        coefficients, IBD sharing probabilities, and inbreeding coefficients. 
        Additionally, functions are provided to perform efficient variance 
        component estimation and mixed model association testing for both 
        quantitative and binary phenotypes.

## Bioconductor

### Release version

[http://www.bioconductor.org/packages/release/bioc/html/GENESIS.html](http://www.bioconductor.org/packages/release/bioc/html/GENESIS.html)

* Tutorials
  * [PC-AiR/PC-Relate](http://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/pcair.html)
  * [GWAS](http://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/assoc_test.html)
  * [Sequence data](http://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/assoc_test_seq.html)
* [Reference manual](http://www.bioconductor.org/packages/release/bioc/manuals/GENESIS/man/GENESIS.pdf)
* [News](http://bioconductor.org/packages/release/bioc/news/GENESIS/NEWS)

### Development version

[http://www.bioconductor.org/packages/devel/bioc/html/GENESIS.html](http://www.bioconductor.org/packages/devel/bioc/html/GENESIS.html)

* Tutorials
  * [PC-AiR/PC-Relate](http://bioconductor.org/packages/devel/bioc/vignettes/GENESIS/inst/doc/pcair.html)
  * [GWAS](http://bioconductor.org/packages/devel/bioc/vignettes/GENESIS/inst/doc/assoc_test.html)
  * [Sequence data](http://bioconductor.org/packages/devel/bioc/vignettes/GENESIS/inst/doc/assoc_test_seq.html)
* [Reference manual](http://www.bioconductor.org/packages/devel/bioc/manuals/GENESIS/man/GENESIS.pdf)
* [News](http://bioconductor.org/packages/devel/bioc/news/GENESIS/NEWS)

## Citation

Gogarten, S.M., Sofer, T., Chen, H., Yu, C., Brody, J.A., Thornton, T.A., Rice, K.M., and Conomos, M.P. (2019). Genetic association testing using the GENESIS R/Bioconductor package. *Bioinformatics*. [doi:10.1093/bioinformatics/btz567](http://doi.org/10.1093/bioinformatics/btz567).

## Installation

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("GENESIS")
```
