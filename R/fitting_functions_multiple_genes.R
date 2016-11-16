#' Fit models.
#' 
#' Model fitting functions for \code{regsplice} package.
#' 
#' There are three model fitting functions:
#' 
#' \code{fit_models_reg} fits regularized (lasso) models containing an optimal subset of 
#' exon:condition interaction terms for each gene. The model fitting procedure penalizes 
#' the interaction terms only, so that the main effect terms for exons and samples are 
#' always included. This ensures that the null model is nested, allowing likelihood ratio
#' tests to be calculated.
#' 
#' \code{fit_models_null} fits the null models, which do not contain any interaction
#' terms.
#' 
#' \code{fit_models_GLM} fits full GLMs, which contain all exon:condition interaction
#' terms for each gene.
#' 
#' See \code{\link{create_design_matrix}} for more details about the terms in each model.
#' 
#' The fitting functions fit models for all genes in the data set. The functions are 
#' parallelized using \code{BiocParallel::bplapply} for faster runtime. For 
#' \code{fit_models_reg}, the default number of processor cores is 8, or the maximum 
#' available if less than 8. For \code{fit_models_null} and \code{fit_models_GLM}, the
#' default is one core, since these functions are already extremely fast for most data
#' sets.
#' 
#' A random seed can be provided with the \code{seed} argument, to generate reproducible
#' results.
#' 
#' 
#' @param Y RNA-seq read counts for multiple genes (list of data frames or matrices). 
#'   Created using \code{\link{prepare_data}}, \code{\link{filter_exons}}, and/or
#'   \code{\link{voom_weights}}.
#' @param condition Experimental conditions for each sample (character or numeric vector,
#'   or factor).
#' @param weights Optional exon-level precision weights (list of data frames or 
#'   matrices). These will usually be generated with \code{\link{voom_weights}}, but may 
#'   also be calculated externally.
#' @param alpha Elastic net parameter \code{alpha} for \code{glmnet} model fitting 
#'   functions. Must be between 0 (ridge regression) and 1 (lasso). Default is 1 (lasso).
#'   See \code{glmnet} documentation for more details.
#' @param lambda_choice Parameter to select which optimal lambda value to choose from the
#'   \code{cv.glmnet} cross validation fit. Choices are "lambda.min" (model with minimum 
#'   cross-validated error) and "lambda.1se" (most regularized model with cross-validated
#'   error within one standard error of minimum). Default is "lambda.min". See 
#'   \code{glmnet} documentation for more details.
#' @param n_cores Number of cores for parallel evaluation. For \code{fit_models_reg}, the
#'   default is 8, or the maximum available if less than 8. For \code{fit_models_GLM} and
#'   \code{fit_models_null}, the default is 1, since these functions are already very
#'   fast.
#' @param seed Random seed (integer). Default is NULL. Provide an integer value to set
#'   the random seed for reproducible results.
#' @param progress_bar Whether to display progress bar (\code{fit_models_reg} only). 
#'   Default is TRUE.
#' @param return_fitted Whether to return fitted model objects. Default is FALSE.
#' @param ... Other arguments to pass to \code{cv.glmnet}, \code{glmnet}, or \code{glm}.
#'   
#' @return Returns a list containing:
#' \itemize{
#' \item dev: Deviance of fitted models for each gene.
#' \item df: Degrees of freedom of fitted models for each gene.
#' \item fit: Fitted model objects (if \code{return_fitted = TRUE}).
#' }
#' 
#' @family create_design_matrix fit_models_reg fit_models_null fit_models_GLM LR_tests
#' 
#' @seealso \code{\link[glmnet]{glmnet}} \code{\link[glmnet]{cv.glmnet}} 
#'   \code{\link[stats]{glm}}
#'   
#' @importFrom BiocParallel multicoreWorkers MulticoreParam bplapply
#' 
#' @export
#' 
#' @examples
#' counts <- matrix(sample(100:200, 40 * 6, replace = TRUE), nrow = 40)
#' gene <- rep(paste0("gene", 1:4), times = c(11, 6, 8, 15))
#' condition <- rep(c(0, 1), each = 3)
#' 
#' Y <- prepare_data(counts, gene)
#' Y <- filter_exons(Y)
#' 
#' # optional 'voom' weights
#' out_voom <- voom_weights(Y, condition)
#' weights <- out_voom$weights
#' 
#' fit_reg  <- fit_models_reg(Y, condition, weights, n_cores = 1)
#' fit_null <- fit_models_null(Y, condition, weights)
#' fit_GLM  <- fit_models_GLM(Y, condition, weights)
#' 
fit_models_reg <- function(Y, condition, weights = NULL, alpha = 1, 
                           lambda_choice = c("lambda.min", "lambda.1se"), 
                           n_cores = NULL, seed = NULL, progress_bar = TRUE, 
                           return_fitted = FALSE, ...) {
  
  if (!(is.list(Y) & !is.data.frame(Y))) {
    stop("data Y for multiple genes must be a list of data frames or matrices")
  }
  
  lambda_choice <- match.arg(lambda_choice)
  
  FUN <- function(i) {
    fit_reg_single(Y = Y[[i]], condition = condition, weights = weights[[i]], 
                   alpha = alpha, lambda_choice = lambda_choice, ...)
  }
  
  message("Fitting regularized (lasso) models...")
  if (is.null(n_cores)) n_cores <- min(BiocParallel::multicoreWorkers(), 8)
  BPPARAM <- BiocParallel::MulticoreParam(workers = n_cores, 
                                          RNGseed = seed, 
                                          progressbar = progress_bar)
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



#' @rdname fit_models_reg
#' @export
#' 
fit_models_null <- function(Y, condition, weights = NULL, 
                            n_cores = 1, seed = NULL, return_fitted = FALSE, ...) {
  
  if (!(is.list(Y) & !is.data.frame(Y))) {
    stop("data Y for multiple genes must be a list of data frames or matrices")
  }
  
  FUN <- function(i) {
    fit_null_single(Y = Y[[i]], condition = condition, weights = weights[[i]], ...)
  }
  
  message("Fitting null models...")
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



#' @rdname fit_models_reg
#' @export
#' 
fit_models_GLM <- function(Y, condition, weights = NULL, 
                           n_cores = 1, seed = NULL, return_fitted = FALSE, ...) {
  
  if (!(is.list(Y) & !is.data.frame(Y))) {
    stop("data Y for multiple genes must be a list of data frames or matrices")
  }
  
  FUN <- function(i) {
    fit_GLM_single(Y = Y[[i]], condition = condition, weights = weights[[i]], ...)
  }
  
  message("Fitting full GLMs...")
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


