#####################################
## Functions for regsplice package ##
## author: Lukas Weber             ##
#####################################


#' Fit multiple genes (parallelized).
#' 
#' Internal function to apply model fitting functions to multiple genes using
#' parallelization.
#' 
#' The internal model fitting functions (\code{\link{fitRegModelSingle}}, 
#' \code{\link{fitGLMSingle}}, and \code{\link{fitNullModelSingle}}) are designed for 
#' fitting single genes only. \code{fitParallel} applies these functions to multiple
#' genes, using parallelization to speed up runtime.
#' 
#' In most cases, the user will access only the external functions
#' \code{\link{fitRegModel}}, \code{\link{fitGLM}}, and \code{\link{fitNull}}. See
#' documentation for those functions for more details.
#' 
#' @param fitting_function Which single-gene fitting function to use. Choices are 
#'   "regularized" (\code{\link{fitRegModelSingle}}), "GLM" (\code{\link{fitGLMSingle}}),
#'   or "null" (\code{\link{fitNullSingle}}). All other arguments (except
#'   \code{fitting_function} and \code{ncores}) will be passed to the single-gene fitting
#'   function. See documentation for \code{\link{fitRegModel}}, \code{\link{fitGLM}}, and
#'   \code{\link{fitNull}} for information about the other arguments.
#' @param ncores Number of cores for parallel evaluation. Default is 1. Note that on
#'   Windows systems, the code will always run on a single core.
#' 
#' @return Returns a list containing:
#' \itemize{
#' \item fit_genes: (list) fitted model objects
#' \item dev_genes: (vector) deviance of fitted models
#' \item df_genes: (vector) degrees of freedom of fitted models
#' }
#' 
#' @seealso \code{\link{fitRegModel}} \code{\link{fitGLM}} \code{\link{fitNullModel}} 
#'   \code{\link{fitRegModelSingle}} \code{\link{fitGLMSingle}}
#'   \code{\link{fitNullModelSingle}}
#'   
fitParallel <- function(fitting_function = c("regularized", "GLM", "null"), 
                        X_genes = NULL, Y_genes, weights_genes = NULL, 
                        group = NULL, nexons_genes = NULL, 
                        alpha = 1, lambda_choice = c("lambda.min", "lambda.1se"), 
                        return_fitted = FALSE, 
                        ncores = 1, ...) {
  
  fitting_function <- match.arg(fitting_function)
  
  if (fitting_function == "regularized") {
    FUN <- function(i) {
      fitRegModelSingle(X = X_genes[[i]], Y = Y_genes[[i]], weights = weights_genes[[i]], 
                        group = group, nexons = nexons_genes[i], 
                        alpha = alpha, lambda_choice = lambda_choice, ...)
    }
  } else if (fitting_function == "GLM") {
    FUN <- function(i) {
      fitGLMSingle(X = X_genes[[i]], Y = Y_genes[[i]], weights = weights_genes[[i]], 
                   group = group, nexons = nexons_genes[i], ...)
    }
  } else if (fitting_function == "null") {
    FUN <- function(i) {
      fitNullModelSingle(X = X_genes[[i]], Y = Y_genes[[i]], weights = weights_genes[[i]], 
                         group = group, nexons = nexons_genes[i], ...)
    }
  }
  
  if (!is.null(X_genes) && !is.list(X_genes)) {
    X_genes <- list(X_genes = X_genes)
  }
  if (!is.list(Y_genes)) {
    Y_genes <- list(Y_genes = Y_genes)
  }
  if (!is.null(weights_genes) && !is.list(weights_genes)) {
    weights_genes <- list(weights_genes = weights_genes)
  }
  
  BPPARAM <- BiocParallel::MulticoreParam(workers = ncores)
  n <- length(Y_genes)
  res <- BiocParallel::bplapply(seq_len(n), FUN = FUN, BPPARAM = BPPARAM)
  
  if (return_fitted == TRUE) fit_genes <- sapply(res, "[[", "fit")
  dev_genes <- sapply(res, "[[", "dev")
  df_genes <- sapply(res, "[[", "df")
  
  if (return_fitted == TRUE) {
    return(list(fit_genes = fit_genes, dev_genes = dev_genes, df_genes = df_genes))
  } else {
    return(list(dev_genes = dev_genes, df_genes = df_genes))
  }
}


