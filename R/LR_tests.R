#' Calculate likelihood ratio tests.
#' 
#' Calculate likelihood ratio tests between fitted models and null models.
#' 
#' The regularized (lasso) fitted models contain an optimal subset of exon:condition 
#' interaction terms; the GLM fitted models contain all exon:condition interaction terms;
#' and the null models do not contain any interaction terms. The null models are
#' therefore nested within the fitted models.
#' 
#' The likelihood ratio (LR) tests compare the fitted models against the nested null
#' models.
#' 
#' If the regularized (lasso) model contains at least one exon:condition interaction 
#' term, the LR test compares the lasso model against the null model. However, if the 
#' lasso model contains zero interaction terms, then the lasso and null models are 
#' identical, so the LR test cannot be calculated. The \code{when_null_selected} argument
#' lets the user choose what to do in these cases: either set p-values equal to 1 
#' (\code{when_null_selected = "ones"}); or calculate a LR test using the full GLM
#' containing all exon:condition interaction terms (\code{when_null_selected = "GLM"}),
#' which reduces power due to the larger number of terms, but allows the evidence for
#' differential exon usage among these genes to be distinguished. You can also return
#' \code{NA}s for these genes (\code{when_null_selected = "NA"}).
#' 
#' The default option is \code{when_null_selected = "ones"}. This simply calls all these
#' genes non-significant, which in most cases is sufficient since we are more interested
#' in genes with strong evidence for differential exon usage. However, if it is important
#' to rank the low-evidence genes in your data set, then use the \code{when_null_selected
#' = GLM} option.
#' 
#' If \code{when_null_selected = "ones"}, the full GLM fitted models are not required, so
#' you can set \code{fitted_models_GLM = NULL} (the default).
#' 
#' 
#' @param fitted_models_reg Regularized (lasso) fitted models (output from
#'   \code{\link{fit_reg}}.
#' @param fitted_models_GLM Full GLM fitted models (output from \code{\link{fit_GLM}}).
#'   Not required if \code{when_null_selected = "ones" or "NA"}. Default is \code{NULL}.
#' @param fitted_models_null Null models (output from \code{\link{fit_null}}).
#' @param when_null_selected Which option to use for genes where the lasso model selects 
#'   zero interaction terms, i.e. identical to the null model. Options are \code{"ones"},
#'   \code{"GLM"}, and \code{"NA"}. Default is \code{"ones"}. See below for details.
#' 
#' @return Returns a list containing:
#' \itemize{
#' \item p_vals: p-values
#' \item p_adj: p-values adjusted for multiple testing (Benjamini-Hochberg false
#'   discovery rates)
#' \item LR_stats: likelihood ratio statistics
#' \item df_tests: degrees of freedom of the likelihood ratio tests
#' }
#' 
#' @family create_design_matrix fit_reg fit_GLM fit_null LR_tests
#' 
#' @export
#' 
#' @examples
#' condition <- rep(c(0, 1), each = 3)
#' n_exons <- 10
#' Y <- list(as.data.frame(matrix(sample(100:200, 60, replace = TRUE), nrow = 10)))
#' fitted_models_reg <- fit_reg(Y, condition)
#' fitted_models_GLM <- fit_GLM(Y, condition)
#' fitted_models_null <- fit_null(Y, condition)
#' 
#' LR_tests(fitted_models_reg = fitted_models_reg, 
#'          fitted_models_GLM = NULL, 
#'          fitted_models_null = fitted_models_null, 
#'          when_null_selected = "ones")
#' 
#' LR_tests(fitted_models_reg = fitted_models_reg, 
#'          fitted_models_GLM = fitted_models_GLM, 
#'          fitted_models_null = fitted_models_null, 
#'          when_null_selected = "GLM")
#' 
LR_tests <- function(fitted_models_reg, fitted_models_GLM = NULL, fitted_models_null, 
                     when_null_selected = c("ones", "GLM", "NA")) {
  
  when_null_selected <- match.arg(when_null_selected)
  
  if (is.null(fitted_models_GLM) & when_null_selected == "GLM") {
    stop('fitted_models_GLM must be provided if when_null_selected = "GLM"')
  }
  
  LR_stats <- abs(unlist(fitted_models_reg$dev) - unlist(fitted_models_null$dev))
  df_tests <- abs(unlist(fitted_models_reg$df) - unlist(fitted_models_null$df))
  
  # genes where lasso selected zero interaction terms (equivalent to null model)
  ix_replace <- df_tests == 0
  
  LR_stats <- LR_stats[!ix_replace]
  df_tests <- df_tests[!ix_replace]
  
  p_vals_keep <- pchisq(LR_stats, df_tests, lower.tail=FALSE)
  
  p_vals <- p_adj <- rep(NA, length(fitted_models_reg$dev))
  
  if (when_null_selected == "ones") {
    p_vals[!ix_replace] <- p_vals_keep
    p_vals[ix_replace] <- 1
    
    p_adj[!ix_replace] <- p.adjust(p_vals_keep, method = "fdr")
    p_adj[ix_replace] <- 1
  
  } else if (when_null_selected == "GLM") {
    LR_stats_GLM <- abs(unlist(fitted_models_GLM$dev) - unlist(fitted_models_null$dev))
    df_tests_GLM <- abs(unlist(fitted_models_GLM$df) - unlist(fitted_models_null$df))
    
    p_vals_GLM <- pchisq(LR_stats_GLM, df_tests_GLM, lower.tail=FALSE)
    
    p_vals[!ix_replace] <- p_vals_keep
    p_vals[ix_replace] <- p_vals_GLM[ix_replace]
    
    p_adj <- p.adjust(p_vals, method = "fdr")
    
  } else if (when_null_selected == "NA") {
    p_vals[!ix_replace] <- p_vals_keep
    
    p_adj[!ix_replace] <- p.adjust(p_vals_keep, method = "fdr")
  }
  
  list(p_vals = p_vals, p_adj = p_adj, LR_stats = LR_stats, df_tests = df_tests)
}


