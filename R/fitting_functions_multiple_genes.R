#' Fit models.
#' 
#' Model fitting functions for \emph{regsplice} package.
#' 
#' There are three model fitting functions:
#' 
#' \code{fit_reg} fits regularized (lasso) models containing an optimal subset of 
#' exon:condition interaction terms for each gene. The model fitting procedure penalizes 
#' the interaction terms only, so that the main effect terms for exons and samples are 
#' always included. This ensures that the null model is nested, allowing likelihood ratio
#' tests to be calculated.
#' 
#' \code{fit_GLM} fits full GLMs containing interaction terms for every exon in each
#' gene.
#' 
#' \code{fit_null} fits the null models, which contain zero zero interaction terms.
#' 
#' See \code{\link{create_design_matrix}} for more details about the terms in each model.
#' 
#' The fitting functions fit models for all genes in the data set, and are parallelized 
#' for faster runtime. For \code{fit_reg}, the default number of cores is 8, or the 
#' maximum available if less than 8. For \code{fit_GLM} and \code{fit_null}, the default 
#' is one core, since these functions are already extremely fast; if they take longer 
#' than a few seconds for your data set, it may be beneficial to try increasing the 
#' number of cores. The parallelization is implemented using
#' \code{BiocParallel::bplapply}.
#' 
#' A random number generation seed can be provided with the \code{seed} argument, to
#' generate reproducible results.
#' 
#' 
#' @param Y RNA-seq counts for multiple genes (list of data frames or matrices). Created 
#'   using the data preparation functions \code{\link{split_genes}} and 
#'   \code{\link{filter_genes}}.
#' @param condition Biological conditions for each sample (character or numeric vector, 
#'   or factor).
#' @param weights Optional weights (list of data frames or matrices), for example 
#'   generated using \code{voom}.
#' @param alpha Elastic net parameter \code{alpha} for \code{glmnet} model fitting 
#'   functions. Must be between 0 (ridge regression) and 1 (lasso). Default is 1 (lasso).
#'   See \code{glmnet} documentation for details.
#' @param lambda_choice Parameter to select which optimal lambda value to choose from the
#'   \code{cv.glmnet} cross validation fit. Choices are "lambda.min" (model with minimum 
#'   cross-validated error) and "lambda.1se" (most regularized model with cross-validated
#'   error within one standard error of minimum). Default is "lambda.min". See
#'   \code{glmnet} documentation for more details.
#' @param return_fitted Whether to return fitted model objects. Default is FALSE.
#' @param n_cores Number of cores for parallel evaluation. For \code{fit_reg}, the 
#'   default is 8, or the maximum available if less than 8. For \code{fit_GLM} and
#'   \code{fit_null}, the default is 1, since these functions are already very fast.
#' @param seed Random number generation seed (integer). Default is NULL. Provide an
#'   integer value to set the random seed for reproducible results.
#' @param ... Other arguments to pass to \code{cv.glmnet}, \code{glmnet}, or \code{glm}.
#'   
#' @return Returns a list containing:
#' \itemize{
#' \item dev Deviance of fitted models for each gene.
#' \item df Degrees of freedom of fitted models for each gene.
#' \item fit Fitted model objects (if \code{return_fitted = TRUE}).
#' }
#' 
#' @family create_design_matrix fit_reg fit_GLM fit_null LR_tests
#' 
#' @seealso \code{\link[glmnet]{glmnet}} \code{\link[glmnet]{cv.glmnet}}
#'   \code{\link[stats]{glm}}
#'   
#' @export
#' 
#' @examples
#' condition <- rep(c(0, 1), each = 3)
#' n_exons <- 10
#' Y <- list(as.data.frame(matrix(sample(100:200, 60, replace = TRUE), nrow = 10)))
#' fit_reg(Y, condition)
#' fit_GLM(Y, condition)
#' fit_null(Y, condition)
#' 
fit_reg <- function(Y, condition, weights = NULL, alpha = 1, 
                    lambda_choice = c("lambda.min", "lambda.1se"), 
                    return_fitted = FALSE, n_cores = NULL, seed = NULL, ...) {
  
  if (!(is.list(Y) & !is.data.frame(Y))) {
    stop("data Y for multiple genes must be a list of data frames or matrices")
  }
  
  lambda_choice <- match.arg(lambda_choice)
  
  FUN <- function(i) {
    fit_reg_single(Y = Y[[i]], condition = condition, weights = weights[[i]], 
                   alpha = alpha, lambda_choice = lambda_choice, ...)
  }
  
  if (is.null(n_cores)) n_cores <- max(BiocParallel::multicoreWorkers(), 8)
  BPPARAM <- BiocParallel::MulticoreParam(workers = n_cores, RNGseed = seed)
  n_genes <- length(Y)
  
  res <- BiocParallel::bplapply(seq_len(n_genes), FUN = FUN, BPPARAM = BPPARAM)
  
  if (return_fitted) fit_genes <- sapply(res, "[[", "fit")
  dev_genes <- sapply(res, "[[", "dev")
  df_genes <- sapply(res, "[[", "df")
  
  if (return_fitted) {
    return(list(dev = dev_genes, df = df_genes, fit = fit_genes))
  } else {
    return(list(dev = dev_genes, df = df_genes))
  }
}



#' @rdname fit_reg
#' @export
#' 
fit_GLM <- function(Y, condition, weights = NULL, 
                    return_fitted = FALSE, n_cores = 1, seed = NULL, ...) {
  
  if (!(is.list(Y) & !is.data.frame(Y))) {
    stop("data Y for multiple genes must be a list of data frames or matrices")
  }
  
  FUN <- function(i) {
    fit_GLM_single(Y = Y[[i]], condition = condition, weights = weights[[i]], ...)
  }
  
  BPPARAM <- BiocParallel::MulticoreParam(workers = n_cores, RNGseed = seed)
  n_genes <- length(Y)
  
  res <- BiocParallel::bplapply(seq_len(n_genes), FUN = FUN, BPPARAM = BPPARAM)
  
  if (return_fitted) fit_genes <- sapply(res, "[[", "fit")
  dev_genes <- sapply(res, "[[", "dev")
  df_genes <- sapply(res, "[[", "df")
  
  if (return_fitted) {
    return(list(dev = dev_genes, df = df_genes, fit = fit_genes))
  } else {
    return(list(dev = dev_genes, df = df_genes))
  }
}



#' @rdname fit_reg
#' @export
#' 
fit_null <- function(Y, condition, weights = NULL, 
                     return_fitted = FALSE, n_cores = 1, seed = NULL, ...) {
  
  if (!(is.list(Y) & !is.data.frame(Y))) {
    stop("data Y for multiple genes must be a list of data frames or matrices")
  }
  
  FUN <- function(i) {
    fit_null_single(Y = Y[[i]], condition = condition, weights = weights[[i]], ...)
  }
  
  BPPARAM <- BiocParallel::MulticoreParam(workers = n_cores, RNGseed = seed)
  n_genes <- length(Y)
  
  res <- BiocParallel::bplapply(seq_len(n_genes), FUN = FUN, BPPARAM = BPPARAM)
  
  if (return_fitted) fit_genes <- sapply(res, "[[", "fit")
  dev_genes <- sapply(res, "[[", "dev")
  df_genes <- sapply(res, "[[", "df")
  
  if (return_fitted) {
    return(list(dev = dev_genes, df = df_genes, fit_genes = fit_genes))
  } else {
    return(list(dev = dev_genes, df = df_genes))
  }
}


