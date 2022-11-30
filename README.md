# ChromENVEE: Chromatin ENVironment and Enhancer-dependent Expression

Standard analyses on ChIPseq data provide information (annotation, enrichment level) at the gene body level but do not necessarily investigate other genomic regions. ChromHMM R package allows to go further by predicting chromatin states using ChIPSeq datasets for several histone marks. The present R package ChromENVEE uses the chromatin states obtained by ChromHMM and compare them with transcriptomic data (RNAseq) and other ChIP-Seq data.

Specifically, ChromENVEE implements functions to associate all the neighbouring genes to a list of enhancers and to define the chromatin environment of genes using chromatin states informations. Several visualization functions are available to summarize the distribution of chromatin states, characterize genes associated with enhancers and also assign chromatin environment to genes.

### Installation

To install ChromENVEE package from GitHub using [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) package:

```
devtools::install_github("ManonCoulee/ChromENVEE")
```
```
library(ChromENVEE)
```

ChromENVEE is available for R  version >= 3.6

### Installation with vignette building

A [vignette](https://github.com/ManonCoulee/ChromENVEE/blob/master/doc/ChromENVEE.pdf) is available to use ChromENVEE.


It is possible to build this vignette in R with the following command:

```
devtools::install_github("ManonCoulee/ChromENVEE", build_vignettes = TRUE)
```

### Citation

To cite the package, please use this following citation:

> Manon Coulee, Guillaume Meurice, Julie Cocquet* and Laila El Khattabi* (2022). ChromENVEE:
Chromatin Environment and Enhancer-dependent Expression. R package version 0.99.8. *co-authorship
