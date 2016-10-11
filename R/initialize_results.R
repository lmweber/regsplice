#' @include class_RegspliceData.R class_RegspliceResults.R
NULL




#' Initialize RegspliceResults object.
#' 
#' Initialize a \code{RegspliceResults} object, which will contain the results of the
#' \code{regsplice} analysis.
#' 
#' Creates a \code{\linkS4class{RegspliceResults}} object containing gene names only. 
#' This object will subsequently be populated using the functions 
#' \code{\link{fit_reg_multiple}}, \code{\link{fit_null_multiple}}, 
#' \code{\link{fit_full_multiple}}, and \code{\link{LR_tests}}.
#' 
#' Previous step: Calculate \code{limma-voom} transformation and weights with 
#' \code{\link{run_voom}}.
#' Next step: Fit models with \code{\link{fit_reg_multiple}},
#' \code{\link{fit_null_multiple}}, and \code{\link{fit_full_multiple}}.
#' 
#' 
#' @param rs_data \code{\linkS4class{RegspliceData}} object. This should contain gene IDs
#'   in a column named \code{gene_IDs} in the row meta-data, which can be accessed with
#'   the accessor function \code{\link{rowData}}.
#' 
#' 
#' @return Returns a \code{\linkS4class{RegspliceResults}} object containing gene IDs
#'   only.
#' 
#' @seealso \code{\linkS4class{RegspliceData}} \code{\linkS4class{RegspliceResults}} 
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
#' rs_data <- filter_zeros(rs_data)
#' rs_data <- filter_low_counts(rs_data)
#' rs_data <- run_normalization(rs_data)
#' rs_data <- run_voom(rs_data)
#' 
#' rs_results <- initialize_results(rs_data)
#' 
initialize_results <- function(rs_data) {
  
  # unique gene IDs in alphabetical order
  gene_IDs <- names(table(rowData(rs_data)$gene_IDs))
  
  RegspliceResults(gene_IDs)
  
}



