#' Calculate likelihood ratio tests.
#' 
#' Calculate likelihood ratio tests between fitted models and null models.
#' 
#' The regularized (lasso) fitted models contain an optimal subset of exon:condition 
#' interaction terms for each gene, and the GLM fitted models contain all exon:condition 
#' interaction terms. The null models contain zero interaction terms, so they are nested
#' within the fitted models.
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
#' to rank the low-evidence genes in your data set, use the \code{when_null_selected = 
#' GLM} option.
#' 
#' If \code{when_null_selected = "ones"} or \code{when_null_selected = "NA"}, the full 
#' GLM fitted models are not required, so you can set \code{fit_GLM = NULL} (the
#' default).
#' 
#' 
#' @param fit_reg Regularized (lasso) fitted models (output from 
#'   \code{\link{fit_models_reg}}).
#' @param fit_null Null models (output from \code{\link{fit_models_null}}).
#' @param fit_GLM Full GLM fitted models (output from \code{\link{fit_models_GLM}}). Not
#'   required if \code{when_null_selected = "ones" or "NA"}. Default is \code{NULL}.
#' @param when_null_selected Which option to use for genes where the lasso model selects 
#'   zero interaction terms, i.e. identical to the null model. Options are \code{"ones"},
#'   \code{"GLM"}, and \code{"NA"}. Default is \code{"ones"}. See below for details.
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
#' @family model-fitting-tests
#' @seealso \code{\link{summary_table}}
#' 
#' @importFrom stats pchisq p.adjust
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
#' LR_tests(fit_reg = fit_reg, 
#'          fit_null = fit_null, 
#'          when_null_selected = "ones")
#' 
#' LR_tests(fit_reg = fit_reg, 
#'          fit_null = fit_null, 
#'          fit_GLM = fit_GLM, 
#'          when_null_selected = "GLM")
#' 
LR_tests <- function(fit_reg, fit_null, fit_GLM = NULL, 
                     when_null_selected = c("ones", "GLM", "NA")) {
  
  if (!identical(fit_reg$gene, fit_null$gene)) {
    stop("gene names do not match (fit_reg and fit_null)")
  }
  if (!is.null(fit_GLM) & !identical(fit_reg$gene, fit_GLM$gene)) {
    stop("gene names do not match (fit_reg and fit_GLM)")
  }
  
  gene <- fit_reg$gene
  
  when_null_selected <- match.arg(when_null_selected)
  
  if (is.null(fit_GLM) & when_null_selected == "GLM") {
    stop('fit_GLM must be provided if when_null_selected = "GLM"')
  }
  
  LR_stats <- abs(unlist(fit_reg$dev) - unlist(fit_null$dev))
  df_tests <- abs(unlist(fit_reg$df) - unlist(fit_null$df))
  
  # genes where lasso selected zero interaction terms (equivalent to null model); 
  # or NAs (where glmnet did not complete)
  ix_remove <- df_tests == 0 | is.na(df_tests)
  
  p_vals_keep <- stats::pchisq(LR_stats[!ix_remove], df_tests[!ix_remove], lower.tail=FALSE)
  
  p_vals <- p_adj <- rep(NA, length(fit_reg$dev))
  
  if (when_null_selected == "ones") {
    p_vals[!ix_remove] <- p_vals_keep
    p_vals[ix_remove] <- 1
    
    # multiple testing adjustment for number of calculated p-values (i.e. non-ones)
    p_adj[!ix_remove] <- stats::p.adjust(p_vals_keep, method = "fdr")
    p_adj[ix_remove] <- 1
    
    LR_stats[ix_remove] <- NA
    df_tests[ix_remove] <- NA
    
  } else if (when_null_selected == "GLM") {
    LR_stats_GLM <- abs(unlist(fit_GLM$dev) - unlist(fit_null$dev))
    df_tests_GLM <- abs(unlist(fit_GLM$df) - unlist(fit_null$df))
    
    p_vals_GLM <- stats::pchisq(LR_stats_GLM, df_tests_GLM, lower.tail=FALSE)
    
    p_vals[!ix_remove] <- p_vals_keep
    p_vals[ix_remove] <- p_vals_GLM[ix_remove]
    
    # multiple testing adjustment for number of calculated p-values (i.e. all genes)
    p_adj <- stats::p.adjust(p_vals, method = "fdr")
    
    LR_stats[ix_remove] <- LR_stats_GLM[ix_remove]
    df_tests[ix_remove] <- df_tests_GLM[ix_remove]
    
  } else if (when_null_selected == "NA") {
    p_vals[!ix_remove] <- p_vals_keep
    
    # multiple testing adjustment for number of calculated p-values (i.e. non-NAs)
    p_adj[!ix_remove] <- stats::p.adjust(p_vals_keep, method = "fdr")
    
    LR_stats[ix_remove] <- NA
    df_tests[ix_remove] <- NA
  }
  
  list(gene = gene, p_vals = p_vals, p_adj = p_adj, LR_stats = LR_stats, df_tests = df_tests)
}


