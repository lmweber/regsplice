#' @include class_RegspliceData.R class_RegspliceResults.R
NULL




#' Wrapper function to run regsplice.
#' 
#' Wrapper function to run a \code{regsplice} analysis with a single command.
#' 
#' This wrapper function runs the \code{regsplice} analysis pipeline with a single 
#' command.
#' 
#' The required input format is a \code{RegspliceData} object, which is created with the
#' \code{\linkS4class{RegspliceData}} constructor function.
#' 
#' The wrapper function calls each of the individual functions in the analysis pipeline 
#' in sequence. You can also run the individual functions directly, which provides 
#' additional flexibility and insight into the statistical methodology. See the vignette 
#' for a description of the individual functions and an example workflow.
#' 
#' After running the analysis pipeline, a summary table of the results can be displayed
#' with \code{\link{summary_table}}.
#' 
#' Note that when using exon microarray data, the filtering, normalization, and 
#' \code{voom} steps should be disabled with the respective arguments.
#' 
#' See \code{\linkS4class{RegspliceData}} for details on constructing the input data
#' object; \code{\link{filter_zeros}} and \code{\link{filter_low_counts}} for details
#' about filtering; \code{\link{run_normalization}} and \code{\link{run_voom}} for
#' details about calculation of normalization factors and \code{voom} transformation and
#' weights; \code{\link{create_design_matrix}} for details about the model design
#' matrices; \code{\link{fit_reg_multiple}}, \code{\link{fit_null_multiple}}, or 
#' \code{\link{fit_full_multiple}} for details about the model fitting functions; and 
#' \code{\link{LR_tests}} for details about the likelihood ratio tests.
#' 
#' 
#' @param rs_data \code{RegspliceData} object containing input data. See
#'   \code{\linkS4class{RegspliceData}} for details.
#' @param filter_zeros Whether to filter zero-count exon bins, using
#'   \code{\link{filter_zeros}}. Default is TRUE. Set to FALSE for exon microarray data.
#' @param filter_low_counts Whether to filter low-count exon bins, using
#'   \code{\link{filter_low_counts}}. Default is TRUE. Set to FALSE for exon microarray
#'   data.
#' @param filter_min_per_exon Filtering parameter for low-count exon bins: minimum number
#'   of reads per exon bin, summed across all biological samples. Default is 6. See
#'   \code{\link{filter_low_counts}} for details.
#' @param filter_min_per_sample Filtering parameter for low-count exon bins: minimum
#'   number of reads per biological sample; i.e. for each exon bin, at least one sample
#'   must have this number of reads. Default is 3. See \code{\link{filter_low_counts}}
#'   for details.
#' @param normalize Whether to calculate normalization factors, using 
#'   \code{\link{run_normalization}}. Default is TRUE. If FALSE, non-normalized library 
#'   sizes will be used. Set to FALSE for exon microarray data.
#' @param norm_method Normalization method to use. Options are \code{"TMM"}, 
#'   \code{"RLE"}, \code{"upperquartile"}, and \code{"none"}. Default is \code{"TMM"}. 
#'   See \code{\link{run_normalization}} for details.
#' @param voom Whether to calculate \code{limma-voom} transformation and weights, using 
#'   \code{\link{run_voom}}. Default is TRUE. If FALSE, model fitting functions will use 
#'   the raw input data (not recommended for count data) with exon bins weighted equally.
#'   Set to FALSE for exon microarray data.
#' @param alpha Elastic net parameter \code{alpha} for \code{glmnet} model fitting 
#'   functions. Must be between 0 (ridge regression) and 1 (lasso). Default is 1 (lasso).
#'   See \code{glmnet} documentation for more details.
#' @param lambda_choice Parameter to select which optimal \code{lambda} value to choose 
#'   from the \code{cv.glmnet} cross validation fit. Choices are "lambda.min" (model with
#'   minimum cross-validated error) and "lambda.1se" (most regularized model with 
#'   cross-validated error within one standard error of minimum). Default is 
#'   "lambda.min". See \code{glmnet} documentation for more details.
#' @param when_null_selected Which option to use for genes where the lasso model selects 
#'   zero interaction terms, i.e. identical to the null model. Options are \code{"ones"},
#'   \code{"full"}, and \code{"NA"}. Default is \code{"ones"}. See \code{\link{LR_tests}}
#'   for details.
#' @param n_cores_reg Number of processor cores for fitting regularized models. Default 
#'   is 8, or the maximum available if less than 8.
#' @param n_cores_null Number of processor cores for fitting null models. Default is 1, 
#'   since this function is already very fast.
#' @param n_cores_full Number of processor cores for fitting full models. Default is 1,
#'   since this function is already very fast.
#' @param seed Random seed (integer). Default is NULL. Provide an integer value to set 
#'   the random seed for reproducible results.
#' @param progress_bar Whether to display progress bar during model fitting (regularized 
#'   models only). Default is TRUE.
#' @param ... Other arguments to pass to \code{cv.glmnet}, \code{glmnet}, or \code{glm}.
#' 
#' 
#' @return Returns a \code{\linkS4class{RegspliceResults}} object containing fitted model
#'   objects and likelihood ratio (LR) test results. The LR test results consist of the
#'   following entries for each gene:
#' \itemize{
#' \item p_vals: raw p-values
#' \item p_adj: multiple testing adjusted p-values (Benjamini-Hochberg false discovery 
#' rates, FDR)
#' \item LR_stats: likelihood ratio test statistics
#' \item df_tests: degrees of freedom of likelihood ratio tests
#' }
#' 
#' 
#' @seealso \code{\linkS4class{RegspliceData}} \code{\linkS4class{RegspliceResults}}
#'   \code{\link{initialize_results}} \code{\link{filter_zeros}}
#'   \code{\link{filter_low_counts}} \code{\link{run_normalization}}
#'   \code{\link{run_voom}} \code{\link{create_design_matrix}}
#'   \code{\link{fit_reg_multiple}} \code{\link{fit_null_multiple}}
#'   \code{\link{fit_full_multiple}} \code{\link{LR_tests}} \code{\link{summary_table}}
#' 
#' @export
#'
#' @examples
#' file_counts <- system.file("extdata/vignette_counts.txt", package = "regsplice")
#' data <- read.table(file_counts, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#' head(data)
#' 
#' counts <- data[, 2:7]
#' tbl_exons <- table(sapply(strsplit(data$exon, ":"), function(s) s[[1]]))
#' gene_IDs <- names(tbl_exons)
#' n_exons <- unname(tbl_exons)
#' condition <- rep(c("untreated", "treated"), each = 3)
#' 
#' rs_data <- RegspliceData(counts, gene_IDs, n_exons, condition)
#' 
#' res <- regsplice(rs_data)
#' 
#' summary_table(res)
#' 
regsplice <- function(rs_data, 
                      filter_zeros = TRUE, filter_low_counts = TRUE, 
                      filter_min_per_exon = 6, filter_min_per_sample = 3, 
                      normalize = TRUE, norm_method = "TMM", voom = TRUE, 
                      alpha = 1, lambda_choice = c("lambda.min", "lambda.1se"), 
                      when_null_selected = c("ones", "full", "NA"), 
                      n_cores_reg = NULL, n_cores_null = 1, n_cores_full = 1, 
                      seed = NULL, progress_bar = TRUE, ...) {
  
  lambda_choice <- match.arg(lambda_choice)
  when_null_selected <- match.arg(when_null_selected)
  
  Y <- rs_data
  
  if (filter_zeros) Y <- filter_zeros(data = Y)
  
  if (filter_low_counts) {
    Y <- filter_low_counts(data = Y, 
                           filter_min_per_exon = filter_min_per_exon, 
                           filter_min_per_sample = filter_min_per_sample)
  }
  
  if (normalize) Y <- run_normalization(data = Y, norm_method = norm_method)
  
  if (voom) Y <- run_voom(data = Y)
  
  res <- initialize_results(data = Y)
  
  res <- fit_reg_multiple(results = res, data = Y, 
                          alpha = alpha, lambda_choice = lambda_choice, 
                          n_cores = n_cores_reg, seed = seed, 
                          progress_bar = progress_bar, ...)
  
  res <- fit_null_multiple(results = res, data = Y, 
                           n_cores = n_cores_null, seed = seed, ...)
  
  if (when_null_selected == "full") {
    res <- fit_full_multiple(results = res, data = Y, 
                             n_cores = n_cores_full, seed = seed, ...)
  }
  
  res <- LR_tests(results = res, when_null_selected = when_null_selected)
  
  res
}



