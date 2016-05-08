regsplice
=========

[![Travis-CI Build Status](https://travis-ci.org/lmweber/regsplice.svg?branch=master)](https://travis-ci.org/lmweber/regsplice)

This repository contains the development version of the R package *regsplice*.

The *regsplice* package implements statistical methods for the detection of differential exon usage (differential splicing) in RNA sequencing (RNA-seq) and microarray data sets.

The *regsplice* methods are based on the use of the lasso (L1-regularization) to improve the power of standard generalized linear models, with fast runtimes compared to other leading approaches. The package will be submitted to [Bioconductor](http://bioconductor.org/); and a paper describing the statistical methodology and comparisons to other methods is currently in preparation.


### How to install

The package can be installed from GitHub:

```{r}
library(devtools)
install_github("lmweber/regsplice")
```


### Dependencies

The *regsplice* package depends on *glmnet* (from CRAN) and *BiocParallel* (from [Bioconductor](http://bioconductor.org/)).

*glmnet* will be installed automatically when you install *regsplice*.

*BiocParallel* will not be installed automatically, since it is from Bioconductor instead of CRAN. It can be installed from Bioconductor with the following code:

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("BiocParallel")
```

