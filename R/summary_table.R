#' @include class_RegspliceData.R class_RegspliceResults.R
NULL




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
#' The main \code{regsplice} functions return results in the form of a 
#' \code{\linkS4class{RegspliceResults}} object, which contains slots for gene names, 
#' fitted model object and results, raw p-values, multiple testing adjusted p-values 
#' (Benjamini-Hochberg FDR), likelihood ratio (LR) test statistics, and degrees of 
#' freedom of the LR tests. See \code{\linkS4class{RegspliceResults}} and the main 
#' \code{regsplice} wrapper function \code{\link{regsplice}} for details.
#' 
#' This function generates a summary table of the results. The results are displayed as a
#' data frame of the top \code{n} most highly significant genes, ranked according to 
#' either FDR or raw p-values, up to a specified significance threshold (e.g. FDR < 
#' 0.05).
#' 
#' The argument \code{rank_by} controls whether to rank by FDR or raw p-values. The
#' default is to rank by FDR.
#' 
#' To display results for all genes up to the significance threshold, set the argument 
#' \code{n = Inf}. To display results for all genes in the data set, set both \code{n = 
#' Inf} and \code{threshold = 1}.
#' 
#' Previous step: Run \code{regsplice} pipeline with the \code{\link{regsplice}} wrapper
#' function (or individual functions up to \code{\link{LR_tests}}).
#' 
#' 
#' @param results \code{\linkS4class{RegspliceResults}} object containing results of a
#'   \code{regsplice} analysis, generated with wrapper function \code{\link{regsplice}}
#'   (or individual functions up to \code{\link{LR_tests}}). See
#'   \code{\linkS4class{RegspliceResults}} for details.
#' @param n Number of genes to display in summary table. Default is 20. If the total
#'   number of significant genes up to the significance threshold is less than \code{n},
#'   only the significant genes are shown. Set to \code{Inf} to display all significant
#'   genes; or set both \code{n = Inf} and \code{threshold = 1} to display all genes in
#'   the data set.
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
#' @seealso \code{\linkS4class{RegspliceResults}} \code{\link{regsplice}}
#' 
#' @importFrom utils head
#' @importFrom methods is
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
#' summary_table(res, n = Inf, threshold = 1)
#' 
summary_table <- function(results, n = 20, threshold = 0.05, 
                          rank_by = c("FDR", "p-value", "none")) {
  
  rank_by <- match.arg(rank_by)
  
  if (!("RegspliceResults" %in% is(results))) {
    stop("'results' must be a 'RegspliceResults' object")
  }
  
  if (rank_by == "FDR") {
    ix <- order(results@p_adj)
  } else if (rank_by == "p-value") {
    ix <- order(results@p_vals)
  } else if (rank_by == "none") {
    ix <- seq_along(results@gene_IDs)
  }
  
  res_display <- data.frame(gene_IDs = results@gene_IDs, 
                            p_vals = results@p_vals, 
                            p_adj = results@p_adj, 
                            LR_stats = results@LR_stats, 
                            df_tests = results@df_tests, 
                            stringsAsFactors = FALSE)
  
  res_ordered <- res_display[ix, ]
  row.names(res_ordered) <- NULL
  
  if (rank_by == "FDR") {
    res_sig <- res_ordered[res_ordered$p_adj <= threshold, ]
  } else if (rank_by == "p-value") {
    res_sig <- res_ordered[res_ordered$p_vals <= threshold, ]
  } else if (rank_by == "none") {
    res_sig <- res_ordered
  }
  
  utils::head(res_sig, n = n)
}



