###########################
### S4 class definition ###
###########################


#' @rdname RegspliceResults
#' @export
#' 
setClass("RegspliceResults", 
         slots = list(gene_IDs = "character", 
                      fit_reg_dev = "numeric", fit_reg_df = "numeric", 
                      fit_null_dev = "numeric", fit_null_df = "numeric", 
                      fit_full_dev = "numeric", fit_full_df = "numeric", 
                      p_vals = "numeric", p_adj = "numeric", 
                      LR_stats = "numeric", df_tests = "numeric"))




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
#' \code{\link{fitRegMultiple}}, \code{\link{fitNullMultiple}},
#' \code{\link{fitFullMultiple}}, and \code{\link{LRTests}}.
#' 
#' The function \code{\link{summaryTable}} can be used to display a summary table of the
#' results.
#' 
#' 
#' @param gene_IDs Gene identifiers or names (character vector).
#' @param x \code{RegspliceResults} object (for accessor functions).
#' 
#' 
#' @field gene_IDs Gene identifiers or names (character vector).
#' @field fit_reg_dev Deviance of fitted regularized (lasso) models from 
#'   \code{\link{fitRegMultiple}}.
#' @field fit_reg_df Degrees of freedom of fitted regularized (lasso) models from 
#'   \code{\link{fitRegMultiple}}.
#' @field fit_null_dev Deviance of fitted null models from 
#'   \code{\link{fitNullMultiple}}.
#' @field fit_null_df Degrees of freedom of fitted null models from 
#'   \code{\link{fitNullMultiple}}.
#' @field fit_full_dev Deviance of fitted full models from 
#'   \code{\link{fitFullMultiple}}.
#' @field fit_full_df Degrees of freedom of fitted full models from 
#'   \code{\link{fitFullMultiple}}.
#' @field p_vals Raw p-values (numeric vector).
#' @field p_adj Multiple testing adjusted p-values (Benjamini-Hochberg false discovery 
#'   rates, FDR).
#' @field LR_stats Likelihood ratio test statistics.
#' @field df_tests Degrees of freedom of likelihood ratio tests.
#' 
#' 
#' @section Accessor functions:
#' 
#' \itemize{
#' \item \code{gene_IDs()}: Accesses gene identifiers or names.
#' \item \code{p_vals()}: Accesses raw p-values.
#' \item \code{p_adj()}: Accesses multiple testing adjusted p-values (Benjamini-Hochberg false discovery rates, FDR).
#' \item \code{LR_stats()}: Accesses likelihood ratio test statistics.
#' \item \code{df_tests()}: Accesses degrees of freedom of likelihood ratio tests.
#' }
#' 
#' 
#' @return Returns an empty \code{RegspliceResults} object.
#' 
#' @seealso \code{\link{fitRegMultiple}} \code{\link{fitNullMultiple}}
#'   \code{\link{fitFullMultiple}} \code{\link{LRTests}} \code{\link{summaryTable}}
#' 
#' @importFrom methods new
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




##########################
### Accessor functions ###
##########################


#' @rdname RegspliceResults
#' @export
#' 
setGeneric("gene_IDs", function(x) {
  standardGeneric("gene_IDs")
})


#' @rdname RegspliceResults
#' @export
#' 
setMethod("gene_IDs", "RegspliceResults", function(x) {
  x@gene_IDs
})




#' @rdname RegspliceResults
#' @export
#' 
setGeneric("p_vals", function(x) {
  standardGeneric("p_vals")
})


#' @rdname RegspliceResults
#' @export
#' 
setMethod("p_vals", "RegspliceResults", function(x) {
  x@p_vals
})




#' @rdname RegspliceResults
#' @export
#' 
setGeneric("p_adj", function(x) {
  standardGeneric("p_adj")
})


#' @rdname RegspliceResults
#' @export
#' 
setMethod("p_adj", "RegspliceResults", function(x) {
  x@p_adj
})




#' @rdname RegspliceResults
#' @export
#' 
setGeneric("LR_stats", function(x) {
  standardGeneric("LR_stats")
})


#' @rdname RegspliceResults
#' @export
#' 
setMethod("LR_stats", "RegspliceResults", function(x) {
  x@LR_stats
})




#' @rdname RegspliceResults
#' @export
#' 
setGeneric("df_tests", function(x) {
  standardGeneric("df_tests")
})


#' @rdname RegspliceResults
#' @export
#' 
setMethod("df_tests", "RegspliceResults", function(x) {
  x@df_tests
})



