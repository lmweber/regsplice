#####################################
## Functions for regsplice package ##
## author: Lukas Weber             ##
#####################################


#' Calculate likelihood ratio test.
#' 
#' Vectorized function to calculate likelihood ratio tests between fitted regularized
#' models (or full models with all interaction terms) and null models.
#' 
#' Calculates likelihood ratio tests to compare the fitted regularized models from 
#' \code{\link{fitRegModel}} (or full models from \code{\link{fitGLM}}) against null 
#' models from \code{\link{fitNullModel}}. The null models do not contain any interaction
#' terms, hence they are nested within the regularized and full models.
#' 
#' This function is vectorized to speed up runtime for data sets with large numbers of
#' genes.
#' 
#' The argument \code{when_null_selected} specifies which method to use for genes where 
#' the regularized model fit is equivalent to the null model, i.e. the lasso fit selected
#' zero interaction terms. The options are:
#' 
#' \itemize{
#' \item "ones": (recommended) Set p-values to 1 for these genes. Under this strategy,
#' the interpretation is that the model found no evidence for differential splicing. This
#' is the simplest and most intuitive strategy, however for some data sets it can result
#' in a large set of genes with indistinguishable levels of evidence.
#' 
#' \item "GLM": Re-fits a full GLM with all interaction terms for these genes. This
#' strategy gives up some power for these genes (since the regularization method is no
#' longer used), in order to obtain a way to rank all genes in the data set by their
#' evidence for differential splicing.
#' 
#' \item "NA": Return NAs for the p-values for these genes. Useful primarily for testing 
#' purposes.
#' }
#' 
#' @param fit_reg Fitted regularized model outputs from \code{\link{fitRegModel}}.
#' @param fit_null Fitted null model outputs from \code{\link{fitNullModel}}.
#' @param fit_GLM Optional fitted full GLM outputs from \code{\link{fitGLM}}. Must be
#'   provided if \code{when_null_selected = "GLM"}.
#' @param when_null_selected Method to use when regularized model is equivalent to the
#'   null model, i.e. the lasso fit selected zero interaction terms. Allowed values are
#'   "ones", "GLM", and "NA" (see above for details).
#' 
#' @return Returns a list containing:
#' \itemize{
#' \item lr_stats: likelihood ratio statistics
#' \item df_tests: degrees of freedom of the likelihood ratio tests
#' \item p_vals: p-values
#' \item p_adj: multiple testing adjusted p-values (false discovery rates)
#' }
#' 
#' @seealso \code{\link{fitRegModel}} \code{\link{fitGLM}} \code{\link{fitNullModel}}
#' 
#' @export
#' 
#' @examples
#' set.seed(1)
#' group <- rep(c(0, 1), each = 3)
#' nexons <- 8
#' X <- createDesignMatrix(group, nexons)
#' Y <- rnorm(nrow(X), mean = 2, sd = 1)
#' ix <- c(7, 8) + (8 * rep(0:5, each = 2))
#' Y[ix] <- Y[ix] + 1
#' fit_reg <- fitRegModel(X, Y)
#' fit_null <- fitNullModel(X, Y)
#' lrTest(fit_reg, fit_null)
#' 
LR_tests <- function(fit_reg, fit_GLM = NULL, fit_null, 
                     when_null_selected = c("ones", "GLM", "NA")) {
  
  when_null_selected <- match.arg(when_null_selected)
  
  if (is.null(fit_GLM) & when_null_selected == "GLM") {
    stop('fit_GLM must be provided with when_null_selected = "GLM"')
  }
  
  LR_stats <- abs(unlist(fit_reg$dev_genes) - unlist(fit_null$dev_genes))
  df_tests <- abs(unlist(fit_reg$df_genes) - unlist(fit_null$df_genes))
  
  # genes where lasso selected the null model
  ix_replace <- df_tests == 0
  
  LR_stats <- LR_stats[!ix_replace]
  df_tests <- df_tests[!ix_replace]
  
  p_vals_keep <- pchisq(LR_stats, df_tests, lower.tail=FALSE)
  
  p_vals <- p_adj <- rep(NA, length(fit_reg$dev_genes))
  
  if (when_null_selected == "ones") {
    p_vals[!ix_replace] <- p_vals_keep
    p_vals[ix_replace] <- 1
    
    p_adj[!ix_replace] <- p.adjust(p_vals_keep, method = "fdr")
    p_adj[ix_replace] <- 1
  
  } else if (when_null_selected == "GLM") {
    LR_stats_GLM_all <- abs(unlist(fit_GLM$dev_genes) - unlist(fit_null$dev_genes))
    df_tests_GLM_all <- abs(unlist(fit_GLM$df_genes) - unlist(fit_null$df_genes))
    p_vals_GLM_all <- pchisq(LR_stats_GLM_all, df_tests_GLM_all, lower.tail=FALSE)
    
    p_vals[!ix_replace] <- p_vals_keep
    p_vals[ix_replace] <- p_vals_GLM_all[ix_replace]
    
    p_adj <- p.adjust(p_vals, method = "fdr")
    
  } else if (when_null_selected == "NA") {
    p_vals[!ix_replace] <- p_vals_keep
    
    p_adj[!ix_replace] <- p.adjust(p_vals_keep, method = "fdr")
  }
  
  list(p_vals = p_vals, p_adj = p_adj, LR_stats = LR_stats, df_tests = df_tests)
}


