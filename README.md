# regsplice

[![Travis-CI Build Status](https://travis-ci.org/lmweber/regsplice.svg?branch=master)](https://travis-ci.org/lmweber/regsplice)
[![codecov](https://codecov.io/gh/lmweber/regsplice/branch/master/graph/badge.svg)](https://codecov.io/gh/lmweber/regsplice)


This repository contains the development version of the R package `regsplice`.

The `regsplice` package implements statistical methods for the detection of differential exon usage (differential splicing) in RNA sequencing (RNA-seq) and exon microarray data sets.

The `regsplice` methods are based on the use of the lasso (L1-regularization) to improve the power of standard generalized linear models. A key advantage is that runtimes are fast compared to other leading approaches. A paper describing the statistical methodology and performance comparisons with other methods is currently in preparation.

The package will be submitted to [Bioconductor](http://bioconductor.org/) for inclusion in the next release.


## How to install

The package can be installed from GitHub with:

```{r}
install.packages("devtools")
library(devtools)
install_github("lmweber/regsplice")
```


## Dependencies

The `regsplice` package depends on:

- `glmnet` (from [CRAN](https://cran.r-project.org/))

- `limma`, `edgeR`, `SummarizedExperiment`, and `BiocParallel` (from [Bioconductor](http://bioconductor.org/))

The `glmnet` package will be installed automatically when you install `regsplice`. The other dependencies will not be installed automatically, since they are from Bioconductor instead of CRAN. They can be installed with:

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite(c("limma", "edgeR", "SummarizedExperiment", "BiocParallel"))
```

