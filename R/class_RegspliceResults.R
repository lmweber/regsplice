###########################
### S4 class definition ###
###########################


#' @rdname RegspliceResults
#' @export
#' 
setClass("RegspliceResults", 
         slots = list(gene_IDs = "character", 
                      fit_reg_models = "list", fit_reg_dev = "numeric", fit_reg_df = "numeric", 
                      fit_null_models = "list", fit_null_dev = "numeric", fit_null_df = "numeric", 
                      fit_full_models = "list", fit_full_dev = "numeric", fit_full_df = "numeric", 
                      p_val = "numeric", p_adj = "numeric", 
                      LR_stat = "numeric", df_test = "numeric"))




############################
### Constructor function ###
############################


#' RegspliceResults objects.
#' 
#' \code{RegspliceResults} objects contain the results of a \code{regsplice} analysis.
#' 
#' \code{RegspliceResults} objects are created with the constructor function
#' \code{RegspliceResults()}, which requires the gene IDs as an argument.
#' 
#' Once created, \code{RegspliceResults} objects are then populated using the functions 
#' \code{\link{fit_reg_multiple}}, \code{\link{fit_null_multiple}},
#' \code{\link{fit_full_multiple}}, and \code{\link{LR_tests}}.
#' 
#' The function \code{\link{summary_table}} can be used to display a summary table of the
#' results.
#' 
#' 
#' @param gene_IDs Gene identifiers or names (character vector).
#' 
#' 
#' @field gene_IDs Gene identifiers or names (character vector).
#' @field fit_reg_models Fitted regularized (lasso) model objects from 
#'   \code{\link{fit_reg_multiple}}.
#' @field fit_reg_dev Deviance of fitted regularized (lasso) models from 
#'   \code{\link{fit_reg_multiple}}.
#' @field fit_reg_df Degrees of freedom of fitted regularized (lasso) models from 
#'   \code{\link{fit_reg_multiple}}.
#' @field fit_null_models Fitted null model objects from \code{\link{fit_null_multiple}}.
#' @field fit_null_dev Deviance of fitted null models from 
#'   \code{\link{fit_null_multiple}}.
#' @field fit_null_df Degrees of freedom of fitted null models from 
#'   \code{\link{fit_null_multiple}}.
#' @field fit_full_models Fitted full model objects from \code{\link{fit_full_multiple}}.
#' @field fit_full_dev Deviance of fitted full models from 
#'   \code{\link{fit_full_multiple}}.
#' @field fit_full_df Degrees of freedom of fitted full models from 
#'   \code{\link{fit_full_multiple}}.
#' @field p_val Raw p-values (numeric vector).
#' @field p_adj Multiple testing adjusted p-values (Benjamini-Hochberg false discovery 
#'   rates, FDR).
#' @field LR_stat Likelihood ratio test statistics.
#' @field df_test Degrees of freedom of likelihood ratio tests.
#' 
#' 
#' @return Returns an empty \code{RegspliceResults} object.
#' 
#' @seealso \code{\link{fit_reg_multiple}} \code{\link{fit_null_multiple}}
#'   \code{\link{fit_full_multiple}} \code{\link{LR_tests}} \code{\link{summary_table}}
#' 
#' @export
#' 
#' @examples
#' # initialize RegspliceResults object
#' gene_IDs <- paste0("gene", 1:5)
#' RegspliceResults(gene_IDs)
#' 
RegspliceResults <- function(gene_IDs) {
  
  new("RegspliceResults", gene_IDs = gene_IDs)
  
}



