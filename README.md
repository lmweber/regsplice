# regsplice

[![Travis-CI Build Status](https://travis-ci.org/lmweber/regsplice.svg?branch=master)](https://travis-ci.org/lmweber/regsplice)
[![codecov](https://codecov.io/gh/lmweber/regsplice/branch/master/graph/badge.svg)](https://codecov.io/gh/lmweber/regsplice)


This repository contains the development version of the R package `regsplice`.

The release version is available from [Bioconductor](https://bioconductor.org/packages/regsplice/).

The `regsplice` package implements statistical methods for the detection of differential exon usage (differential splicing) in RNA sequencing (RNA-seq) and exon microarray data sets. The `regsplice` methods are based on the use of the lasso (L1-regularization) to improve the power of standard generalized linear models. A key advantage is that runtimes are fast compared to other leading approaches.

A paper describing the statistical methodology and performance comparisons with other methods is currently in preparation.


## Install release version

The release version can be installed from [Bioconductor](https://bioconductor.org/packages/regsplice/) using the Bioconductor installer. This will also install all required dependencies. This is the recommended option for most users.

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("regsplice")
```


## Install development version

The development version can be installed using the "Devel" version of Bioconductor (see [Bioconductor documentation](http://bioconductor.org/developers/how-to/useDevel/) for details).

Alternatively, you can also install the development version from this GitHub repository using `devtools::install_github`.

```{r}
install.packages("devtools")
library(devtools)
install_github("lmweber/regsplice")
```


## Dependencies

The `regsplice` package depends on:

- `glmnet`, `pbapply` (from [CRAN](https://cran.r-project.org/))

- `limma`, `edgeR`, `SummarizedExperiment` (from [Bioconductor](http://bioconductor.org/))

If you install using the Bioconductor installer, all dependencies will be installed automatically.

If you install from GitHub, the Bioconductor dependencies need to be installed separately.

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite(c("limma", "edgeR", "SummarizedExperiment"))
```

