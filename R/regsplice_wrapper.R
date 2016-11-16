#' Wrapper function to run regsplice.
#' 
#' Wrapper function to run a complete \code{regsplice} analysis with one command.
#' 
#' This wrapper function runs the complete \code{regsplice} pipeline with a single 
#' command. It calls each of the individual functions within the package in sequence. You
#' can also run the individual functions separately, which provides additional 
#' flexibility and insight into the statistical methodology. See the vignette for a 
#' description of the individual functions and an example workflow.
#' 
#' Required inputs are \code{counts} (matrix or data frame of RNA-seq read counts or exon
#' microarray intensities), \code{gene} (vector of gene IDs), and \code{condition} 
#' (vector of experimental conditions).
#' 
#' The gene ID vector \code{gene} must contain one entry for each exon, with repeated 
#' entries for multiple exons within the same gene, so that its length is equal to the 
#' number of rows in \code{counts}; the repeated entries are used to determine gene 
#' lengths. See the vignette for an example showing how to construct the gene ID vector 
#' from a column of gene:exon IDs.
#' 
#' Exon microarray intensities should be log2-transformed, which can either be done 
#' externally (for example with \code{limma-voom}) or with the \code{voom_norm = TRUE} 
#' argument. Note that filtering parameters \code{filter_n1} and \code{filter_n2} need to
#' be adjusted carefully when using exon microarray data; filtering may also be disabled
#' with \code{filter = FALSE}.
#' 
#' See \code{\link{prepare_data}}, \code{\link{filter_exons}}, and
#' \code{\link{voom_weights}} for details about data preparation, filtering low-count
#' exons, and calculation of exon-level precision weights and/or scale normalization;
#' \code{\link{create_design_matrix}} for details about the model design matrices; 
#' \code{\link{fit_models_reg}}, \code{\link{fit_models_null}}, or 
#' \code{\link{fit_models_GLM}} for details about the model fitting functions; and 
#' \code{\link{LR_tests}} for details about the likelihood ratio tests.
#' 
#' 
#' @param counts RNA-seq read counts or exon microarray intensities (matrix or data
#'   frame). Each row is an exon, and each column is a biological sample.
#' @param gene Vector of gene IDs (character vector). Length is equal to the number of 
#'   rows in \code{counts}.
#' @param condition Experimental conditions for each sample (character or numeric vector, 
#'   or factor).
#' @param voom_weights Whether to use \code{limma-voom} exon-level precision weights. 
#'   Default is TRUE. See \code{\link{voom_weights}} for details. If set to FALSE,
#'   weights are not used; i.e. exons are weighted equally.
#' @param voom_norm Whether to use \code{limma-voom} log2-counts per million continuous 
#'   transformation and scale normalization across samples. Default is FALSE. Should be
#'   set to TRUE if using exon microarray intensities that have not already been
#'   log2-transformed. See \code{\link{voom_weights}} for details.
#' @param alpha Elastic net parameter \code{alpha} for \code{glmnet} model fitting 
#'   functions. Must be between 0 (ridge regression) and 1 (lasso). Default is 1 (lasso).
#'   See \code{glmnet} documentation for more details.
#' @param lambda_choice Parameter to select which optimal lambda value to choose from the
#'   \code{cv.glmnet} cross validation fit. Choices are "lambda.min" (model with minimum 
#'   cross-validated error) and "lambda.1se" (most regularized model with cross-validated
#'   error within one standard error of minimum). Default is "lambda.min". See 
#'   \code{glmnet} documentation for more details.
#' @param when_null_selected Which option to use for genes where the lasso model selects 
#'   zero interaction terms, i.e. identical to the null model. Options are \code{"ones"},
#'   \code{"GLM"}, and \code{"NA"}. Default is \code{"ones"}. See \code{\link{LR_tests}} 
#'   for details.
#' @param n_cores_reg Number of processor cores for fitting regularized models. Default
#'   is 8, or the maximum available if less than 8.
#' @param n_cores_null Number of processor cores for fitting null models. Default is 1,
#'   since this function is already very fast.
#' @param n_cores_GLM Number of processor cores for fitting GLMs. Default is 1, since
#'   this function is already very fast.
#' @param seed Random seed (integer). Default is NULL. Provide an integer value to set
#'   the random seed for reproducible results.
#' @param progress_bar Whether to display progress bar during model fitting (regularized 
#'   models only). Default is TRUE.
#' @param return_fitted Whether to return fitted model objects. Default is FALSE.
#' @param filter Whether to filter low-count exons. Default is TRUE. Set to FALSE to
#'   disable filtering (for example, this may be required when using exon microarray
#'   data).
#' @param filter_n1 Parameter for filtering low-count exons: minimum number of reads per
#'   exon, summed across all biological samples. Default is 6. See
#'   \code{\link{filter_exons}} for more details.
#' @param filter_n2 Parameter for filtering low-count exons: minimum number of reads for
#'   a single biological sample per exon, i.e. at least one sample must have this number
#'   of reads. Default is 3. See \code{\link{filter_exons}} for more details.
#' @param ... Other arguments to pass to \code{cv.glmnet}, \code{glmnet}, or \code{glm}.
#' 
#' 
#' @return Returns a list containing:
#' \itemize{
#' \item gene: gene names
#' \item p_vals: raw p-values
#' \item p_adj: multiple testing adjusted p-values (Benjamini-Hochberg false discovery
#' rates, FDR)
#' \item LR_stats: likelihood ratio test statistics
#' \item df_tests: degrees of freedom of likelihood ratio tests
#' }
#' 
#' @seealso \code{\link{prepare_data}} \code{\link{filter_exons}} 
#'   \code{\link{voom_weights}} \code{\link{create_design_matrix}} 
#'   \code{\link{fit_models_reg}} \code{\link{fit_models_null}} 
#'   \code{\link{fit_models_GLM}} \code{\link{LR_tests}} \code{\link{summary_table}}
#' 
#' @export
#'
#' @examples
#' counts <- matrix(sample(100:200, 65 * 6, replace = TRUE), nrow = 65)
#' gene <- rep(paste0("gene", 1:8), times = c(11, 2, 8, 15, 6, 7, 6, 10))
#' condition <- rep(c(0, 1), each = 3)
#' 
#' regsplice(counts, gene, condition)
#' 
regsplice <- function(counts, gene, condition, 
                      voom_weights = TRUE, voom_norm = FALSE, 
                      alpha = 1, lambda_choice = c("lambda.min", "lambda.1se"), 
                      when_null_selected = c("ones", "GLM", "NA"), 
                      n_cores_reg = NULL, n_cores_null = 1, n_cores_GLM = 1, 
                      seed = NULL, progress_bar = TRUE, return_fitted = FALSE, 
                      filter = TRUE, filter_n1 = 6, filter_n2 = 3, 
                      ...) {
  
  lambda_choice <- match.arg(lambda_choice)
  when_null_selected <- match.arg(when_null_selected)
  
  Y <- prepare_data(counts = counts, gene = gene)
  
  if (filter) Y <- filter_exons(Y = Y, n1 = filter_n1, n2 = filter_n2)
  
  if (voom_weights | voom_norm) {
    out_voom <- voom_weights(Y = Y, condition = condition, norm = TRUE)
  }
  
  if (voom_weights) {
    weights <- out_voom$weights
  } else {
    weights <- NULL
  }
  
  if (voom_norm) Y <- out_voom$Y
  
  fit_reg <- fit_models_reg(Y = Y, condition = condition, weights = weights, 
                            alpha = alpha, lambda_choice = lambda_choice, 
                            n_cores = n_cores_reg, seed = seed, 
                            progress_bar = progress_bar, 
                            return_fitted = return_fitted, ...)
  
  fit_null <- fit_models_null(Y = Y, condition = condition, weights = weights, 
                              n_cores = n_cores_null, seed = seed, 
                              return_fitted = return_fitted, ...)
  
  if (when_null_selected == "GLM") {
    fit_GLM <- fit_models_GLM(Y = Y, condition = condition, weights = weights, 
                              n_cores = n_cores_GLM, seed = seed, 
                              return_fitted = return_fitted, ...)
  } else {
    fit_GLM <- NULL
  }
  
  LR_tests(fit_reg = fit_reg, 
           fit_null = fit_null, 
           fit_GLM = fit_GLM, 
           when_null_selected = when_null_selected)
}


