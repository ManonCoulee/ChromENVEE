# ChromENVEE (Chromatin ENVironment and Enhancer Expression)

This package implements functions to associate genes with enhancers, define the chromatin environment of the gene from genomic data (e.g., ChromHMM output or a bed file). Several visualization functions are available to summarize the distribution of chromatin states, characterize genes associated with enhancers and also estimate the chromatin environment of genes.

### Installation

To install ChromENVEE package from GitHub using [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) package :

```
devtools::install_github("ManonCoulee/ChromENVEE")
```
```
library(ChromENVEE)
```

ChromENVEE is available for R  version >= 3.6

### Installation with vignette building

A [vignette](https://github.com/ManonCoulee/ChromENVEE/blob/master/doc/ChromENVEE.pdf) is available to use ChromENVEE
It's possible to build this vignette in R with the following command :

```
devtools::install_github("ManonCoulee/ChromENVEE", build_vignettes = TRUE)
```

### Citation

To cite the package, please use this following citation :
