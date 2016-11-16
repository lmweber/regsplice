#' Summary table.
#' 
#' Display summary table of results from a \code{regsplice} analysis.
#' 
#' The results of a \code{regsplice} analysis consist of a set of multiple testing 
#' adjusted p-values (Benjamini-Hochberg false discovery rates, FDR) quantifying the 
#' statistical evidence for differential exon usage (DEU) for each gene. Typically, the
#' adjusted p-values are used to rank the genes in the data set according to their
#' evidence for DEU, and an appropriate significance threshold (e.g. FDR < 0.05) can be
#' used to generate a list of genes with statistically significant evidence for DEU.
#' 
#' The main \code{regsplice} functions return results structured as a list, containing 
#' gene names, raw p-values, multiple testing adjusted p-values (Benjamini-Hochberg FDR),
#' likelihood ratio (LR) test statistics, and degrees of freedom of the LR tests. See the
#' wrapper function \code{\link{regsplice}} for details.
#' 
#' This function generates a summary table of the results. The results are displayed as a
#' data frame of the top \code{n} most highly significant genes, ranked according to 
#' either FDR or raw p-values, up to a specified significance threshold (e.g. FDR <
#' 0.05).
#' 
#' The argument \code{rank_by} controls whether to rank by FDR or raw p-values.
#' 
#' To display results for all genes up to the significance threshold, set the argument
#' \code{n = Inf}. To display results for all genes in the data set, set both \code{n =
#' Inf} and \code{threshold = 1}.
#' 
#' 
#' @param res Results object containing results of a \code{regsplice} analysis, generated
#'   by either \code{\link{regsplice}} (wrapper function) or \code{\link{LR_tests}}. See 
#'   \code{\link{regsplice}} or \code{\link{LR_tests}} for details.
#' @param n Number of genes to display. Default is 20. If the total number of significant
#'   genes up to the significance threshold is less than \code{n}, only the significant
#'   genes are shown. Set to \code{Inf} to display all significant genes; or set both
#'   \code{n = Inf} and \code{threshold = 1} to display all genes in the data set.
#' @param threshold Significance threshold (for either FDR or raw p-values, depending on 
#'   choice of argument \code{rank_by}). Default is 0.05. Set to 1 to display all
#'   \code{n} genes; or set both \code{n = Inf} and \code{threshold = 1} to display all
#'   genes in the data set.
#' @param rank_by Whether to rank genes by false discovery rate (FDR), raw p-values, or 
#'   no ranking. Choices are \code{"FDR"}, \code{"p-value"}, and \code{"none"}. Default 
#'   is \code{"FDR"}.
#' 
#' 
#' @return Returns a data frame containing results for the top \code{n} most highly 
#'   significant genes, up to the specified significance threshold for the FDR or raw
#'   p-values.
#' 
#' @seealso \code{\link{regsplice}}
#' 
#' @importFrom utils head
#' 
#' @export
#' 
#' @examples
#' file_counts <- system.file("extdata/vignette_counts.txt", package = "regsplice")
#' data <- read.table(file_counts, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#' counts <- data[, 2:7]
#' gene <- sapply(strsplit(data$exon, ":"), function(s) s[[1]])
#' condition <- rep(c("untreated", "treated"), each = 3)
#' 
#' res <- regsplice(counts, gene, condition)
#' 
#' summary_table(res)
#' summary_table(res, n = Inf, threshold = 1)
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
    ix <- order(res$p_val)
  } else if (rank_by == "none") {
    ix <- 1:length(res$p_adj)
  }
  
  ordered <- as.data.frame(res, stringsAsFactors = FALSE)[ix, ]
  row.names(ordered) <- NULL
  
  if (rank_by == "FDR") {
    sig <- ordered[ordered$p_adj <= threshold, ]
  } else if (rank_by == "p-value") {
    sig <- ordered[ordered$p_val <= threshold, ]
  } else if (rank_by == "none") {
    sig <- ordered
  }
  
  utils::head(sig, n = n)
}

