#' Wrapper function to run regsplice.
#' 
#' Wrapper function to run a complete regsplice analysis with one command.
#' 
#' This wrapper function runs the complete \emph{regsplice} pipeline with a single 
#' command. It calls each of the individual functions within the package in sequence. You
#' can also run the individual functions separately, which provides additional 
#' flexibility and insight into the statistical methodology. See the vignette for a 
#' description of the individual functions and an example workflow.
#' 
#' Required inputs are \code{counts} (matrix of RNA-seq counts), \code{gene} (vector of 
#' gene IDs), and \code{condition} (vector of biological conditions).
#' 
#' The gene ID vector \code{gene} must contain repeated entries for multi-exon genes, so 
#' that its length is equal to the total number of rows in \code{counts}; the gene 
#' lengths are taken from the number of repeated entries for each gene. Usually you will 
#' be able to construct the gene ID vector from the row names of a raw data frame. See
#' the vignette for an example.
#' 
#' See \code{\link{create_design_matrix}} for details about the model design matrices; 
#' \code{\link{fit_reg}}, \code{\link{fit_GLM}}, or \code{\link{fit_null}} for details
#' about the model fitting functions; and \code{\link{LR_tests}} for details about the
#' likelihood ratio tests.
#' 
#' 
#' @param counts RNA-seq counts (matrix or data frame). Each row is an exon, and each 
#'   column is a biological sample.
#' @param gene Vector of gene IDs (character vector). Length is equal to the number of
#'   rows in \code{counts}.
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
#' @param when_null_selected Which option to use for genes where the lasso model selects 
#'   zero interaction terms, i.e. identical to the null model. Options are \code{"ones"},
#'   \code{"GLM"}, and \code{"NA"}. Default is \code{"ones"}. See \code{\link{LR_tests}}
#'   for details.
#' @param return_fitted Whether to return fitted model objects. Default is FALSE.
#' @param n_cores_reg Number of cores for fitting regularized models. Default is 8, or
#'   the maximum available if less than 8.
#' @param n_cores_GLM Number of cores for fitting GLMs. Default is 1, since this function
#'   is already very fast.
#' @param n_cores_null Number of cores for fitting null models. Default is 1, since this
#'   function is already very fast.
#' @param seed Random number generation seed (integer). Default is NULL. Provide an
#'   integer value to set the random seed for reproducible results.
#' @param ... Other arguments to pass to \code{cv.glmnet}, \code{glmnet}, or \code{glm}.
#' 
#' 
#' @return Returns a list containing:
#' \itemize{
#' \item dev Deviance of fitted models for each gene.
#' \item df Degrees of freedom of fitted models for each gene.
#' \item fit Fitted model objects (if \code{return_fitted = TRUE}).
#' }
#' 
#' @seealso \code{\link{split_genes}} \code{\link{filter_zeros}}
#'   \code{\link{filter_single_exons}} \code{\link{create_design_matrix}}
#'   \code{\link{fit_reg}} \code{\link{fit_GLM}} \code{\link{fit_null}}
#'   \code{\link{LR_tests}}
#' 
#' @export
#'
#' @examples
#' n_exons <- c(4, 11, 3, 2, 5, 5, 7, 9, 4, 20, 18, 3, 9)
#' n_genes <- length(n_exons)
#' gene <- paste0("gene", rep(1:n_genes, times = n_exons))
#' condition <- rep(c(0, 1), each = 3)
#' 
#' counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
#' 
#' regsplice(counts, gene, condition)
#' 
regsplice <- function(counts, gene, condition, weights = NULL, 
                      alpha = 1, lambda_choice = c("lambda.min", "lambda.1se"), 
                      when_null_selected = c("ones", "GLM", "NA"), 
                      return_fitted = FALSE, 
                      n_cores_reg = NULL, n_cores_GLM = 1, n_cores_null = 1, 
                      seed = NULL, ...) {
  
  lambda_choice <- match.arg(lambda_choice)
  when_null_selected <- match.arg(when_null_selected)
  
  Y <- split_genes(counts = counts, gene = gene)
  Y <- filter_zeros(Y = Y)
  Y <- filter_single_exons(Y = Y)
  
  fitted_models_reg <- fit_reg(Y = Y, condition = condition, weights = weights, 
                               alpha = alpha, lambda_choice = lambda_choice, 
                               return_fitted = return_fitted, 
                               n_cores = n_cores_reg, seed = seed, ...)
  
  if (when_null_selected == "GLM") {
    fitted_models_GLM <- fit_GLM(Y = Y, condition = condition, weights = weights, 
                                 return_fitted = return_fitted, 
                                 n_cores = n_cores_GLM, seed = seed, ...)
  } else {
    fitted_models_GLM <- NULL
  }
  
  fitted_models_null <- fit_null(Y = Y, condition = condition, weights = weights, 
                                 return_fitted = return_fitted, 
                                 n_cores = n_cores_null, seed = seed, ...)
  
  LR_tests(fitted_models_reg = fitted_models_reg, 
           fitted_models_GLM = fitted_models_GLM, 
           fitted_models_null = fitted_models_null, 
           when_null_selected = when_null_selected)
}


