# regsplice

[![Travis-CI Build Status](https://travis-ci.org/lmweber/regsplice.svg?branch=master)](https://travis-ci.org/lmweber/regsplice)

[![codecov](https://codecov.io/gh/lmweber/regsplice/branch/master/graph/badge.svg)](https://codecov.io/gh/lmweber/regsplice)


This repository contains the development version of the R package `regsplice`.

The `regsplice` package implements statistical methods for the detection of differential exon usage (differential splicing) in RNA sequencing (RNA-seq) and exon microarray data sets.

The `regsplice` methods are based on the use of the lasso (L1-regularization) to improve the power of standard generalized linear models. A key advantage is fast runtimes compared to other leading approaches. A paper describing the detailed statistical methodology and performance comparisons with other methods is currently in preparation.

The package will be submitted to [Bioconductor](http://bioconductor.org/) for inclusion in the next release.


## How to install

The package can be installed from GitHub as follows:

```{r}
install.packages("devtools")
library(devtools)
install_github("lmweber/regsplice")
```


## Dependencies

The `regsplice` package depends on `glmnet` (from CRAN) and `BiocParallel` (from [Bioconductor](http://bioconductor.org/)).

`glmnet` will be installed automatically when you install `regsplice`.

`BiocParallel` will not be installed automatically, since it is from Bioconductor instead of CRAN. It can be installed with:

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("BiocParallel")
```

