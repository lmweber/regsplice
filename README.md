regsplice
=========

This repository contains the development version of the R package `regsplice`.

The `regsplice` package implements methods to rank genes in RNA-seq or microarray data sets by their evidence for differential exon usage (splicing), using L1 regularization (lasso) for model selection to improve power.

The package will be submitted to [Bioconductor](http://bioconductor.org/), and a preprint of the accompanying paper will be made available on [bioRxiv](http://biorxiv.org/).


### How to install

The package can be installed from GitHub using the following code. The authentication token is required since this is a private repository.

```{r}
library(devtools)
auth_token <- "ad0c4ab9b6216af8a37b440a26fa2e711270814e"
devtools::install_github("lmweber/regsplice", auth_token = auth_token)
```

### `glmnet` version

(remove this if `glmnet` has been updated)

`regsplice` depends on the R package `glmnet` for the main lasso model fitting procedures. However, there is a bug in the current version of `glmnet` on CRAN (version 2.0-2), which causes errors with the `regsplice` functions. The previous version (1.9-8) works fine, and can be installed from source with the code below.

```{r}
url <- "https://cran.r-project.org/src/contrib/Archive/glmnet/glmnet_1.9-8.tar.gz"
devtools::install_url(url, repos = NULL, type = "source")
```

The `glmnet` authors have informed that the bug will be fixed in the next version of `glmnet`, so in future this step will no longer be required.

See the repository [glmnet-error-example](https://github.com/lmweber/glmnet-error-example) for more details about the bug.


### Dependencies

The `regsplice` package depends on `glmnet` (from CRAN) and `BiocParallel` (from Bioconductor).

(remove this if `glmnet` has been updated:) In the `DESCRIPTION` file, these packages are listed in the `Suggests` field, not `Depends`. This is so that they do not install automatically when you install the `regsplice` package. This is required so that the older version of `glmnet` (1.9-8) can be installed separately from source. As mentioned above, the current version of `glmnet` from CRAN (2.0-2) has a bug, so version 1.9-8 should be used instead.

`BiocParallel` is used to parallelize the `regsplice` functions. The parallelization implemented in `regsplice` works on Mac OS X and Linux systems. On Windows systems, the code will run on a single core (this may change with the upcoming release of Bash and the Linux command line for Windows 10, [expected mid-2016](http://www.theverge.com/2016/3/30/11331014/microsoft-windows-linux-ubuntu-bash)).

`BiocParallel` can be installed from Bioconductor as usual with:

```{r}
source("http://bioconductor.org/biocLite.R")
biocLite("BiocParallel")
```

