# ChromENVEE: Chromatin ENVironment and Enhancer-dependent Expression

ChromENVEE is a package developed to study chromatin states.

This package implements functions to associate all the neighbouring genes to a list of enhancers and to define the chromatin environment of genes using chromatin states informations (e.g., ChromHMM output). Several visualization functions are available to summarize the distribution of chromatin states, characterize genes associated with enhancers and also assign chromatin environment to genes.

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

```
Manon Coulee, Guillaume Meurice, Mitra Barzine, Laila El Khattabi and Julie Cocquet (2022).     ChromENVEE: Chromatin Environment and Enhancer-dependent Expression. R package version 1.1.8.
```
