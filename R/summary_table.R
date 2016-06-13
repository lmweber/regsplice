#' Summary table.
#' 
#' Display summary table of results from a \code{regsplice} analysis.
#' 
#' The results of a \code{regsplice} analysis consist of a set of multiple testing 
#' adjusted p-values (Benjamini-Hochberg false discovery rates, FDR) quantifying the 
#' statistical evidence for differential exon usage (DEU) for each gene. Typically, the 
#' adjusted p-values will be used to rank the genes in the data set according to their 
#' evidence for DEU, and an appropriate significance threshold (e.g. FDR < 0.05) is used 
#' to generate a list of genes with statistically significant evidence for DEU.
#' 
#' The main \code{regsplice} functions return results structured as a list, containing 
#' gene names, raw p-values, multiple testing adjusted p-values (Benjamini-Hochberg FDR),
#' likelihood ratio (LR) test statistics, and degrees of freedom of the LR tests. See the
#' wrapper function \code{\link{regsplice}} for details.
#' 
#' This function converts the results to a more readable format. The results are
#' displayed as a data frame containing the top \code{n} most highly significant genes, 
#' ranked according to either FDR or raw p-values, up to a specified significance
#' threshold.
#' 
#' Set \code{n = Inf} to display results for all genes up to the specified significance 
#' threshold.
#' 
#' Set \code{n = Inf} and \code{threshold = 1} to display results for all genes in the 
#' data set.
#' 
#' 
#' @param res Results object containing results of a \code{regsplice} analysis, generated
#'   by either \code{\link{regsplice}} (wrapper function) or \code{\link{LR_tests}}. See 
#'   \code{\link{regsplice}} or \code{\link{LR_tests}} for details.
#' @param n Number of genes to display. Default is 20. Set to \code{Inf} to display all 
#'   genes up to the specified significance threshold.
#' @param threshold Significance threshold. Default is 0.05. Set to 1 to display all 
#'   genes.
#' @param rank_by Whether to rank genes by false discovery rate (FDR), raw p-values, or 
#'   no ranking. Choices are \code{"FDR"}, \code{"p-value"}, and \code{"none"}. Default 
#'   is \code{"FDR"}.
#' 
#' 
#' @return Returns a data frame containing results for the top \code{n} most highly 
#'   significant genes, up to the specified significance threshold.
#' 
#' @seealso \code{\link{regsplice}} \code{\link{LR_tests}}
#' 
#' @importFrom utils head
#' 
#' @export
#' 
#' @examples
#' counts <- matrix(sample(100:200, 65 * 6, replace = TRUE), nrow = 65)
#' gene <- rep(paste0("gene", 1:8), times = c(11, 2, 8, 15, 6, 7, 6, 10))
#' condition <- rep(c(0, 1), each = 3)
#' 
#' res <- regsplice(counts, gene, condition)
#' 
#' summary_table(res)
#' 
summary_table <- function(res, n = 20, threshold = 0.05, 
                          rank_by = c("FDR", "p-value", "none")) {
  
  if (!is.list(res)) stop("results object must be a list")
  
  if (length(unique(sapply(res, length))) != 1) {
    stop("all list items in results object must have equal length")
  }
  
  rank_by <- match.arg(rank_by)
  
  if (rank_by == "FDR") {
    ix <- order(res$p_adj)
  } else if (rank_by == "p-value") {
    ix <- order(res$p_vals)
  } else if (rank_by == "none") {
    ix <- 1:length(res$p_adj)
  }
  
  ordered <- as.data.frame(res, stringsAsFactors = FALSE)[ix, ]
  row.names(ordered) <- NULL
  
  if (rank_by == "FDR") {
    sig <- ordered[ordered$p_adj <= threshold, ]
  } else if (rank_by == "p-value") {
    sig <- ordered[ordered$p_vals <= threshold, ]
  } else if (rank_by == "none") {
    sig <- ordered
  }
  
  head(sig, n = n)
}

