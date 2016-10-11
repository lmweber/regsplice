#' @include class_RegspliceData.R class_RegspliceResults.R
NULL




#' Calculate 'voom' transformation and weights.
#' 
#' Use \code{limma-voom} to transform counts and calculate exon-level weights.
#' 
#' Raw counts do not fulfill the statistical assumptions required for linear modeling. 
#' The \code{limma-voom} methodology transforms counts to log2-counts per million 
#' (logCPM), and calculates exon-level weights based on the observed mean-variance 
#' relationship. Linear modeling methods can then be applied.
#' 
#' For more details, see the documentation for \code{\link[limma]{voom}} in the 
#' \code{limma} package.
#' 
#' Note that \code{voom} assumes that exon bins (rows) with zero or low counts have 
#' already been removed, so this step should be done after filtering with 
#' \code{\link{filterZeros}} and \code{\link{filterLowCounts}}.
#' 
#' Normalization factors can be provided in a column named \code{norm_factors} in the 
#' column meta-data (\code{colData} slot) of the \code{\linkS4class{RegspliceData}} 
#' object. These will be used by \code{voom} to calculate normalized library sizes. If 
#' normalization factors are not provided, \code{voom} will use non-normalized library 
#' sizes (columnwise total counts) instead.
#' 
#' The experimental conditions or group labels for each biological sample are assumed to 
#' be in a column named \code{condition} in the column meta-data (\code{colData} slot) of
#' the \code{\linkS4class{RegspliceData}} object. This column is created when the object 
#' is initialized with the \code{RegspliceData()} constructor function.
#' 
#' The transformed counts are stored in the updated \code{counts} matrix, which can be 
#' accessed with the \code{\link{countsData}} accessor function. The weights are stored 
#' in a new data matrix labeled \code{weights}, which can be accessed with the 
#' \code{\link{weightsData}} accessor function. In addition, the normalized library sizes
#' (if available) are stored in a new column named \code{lib_sizes} in the column 
#' meta-data (\code{colData} slot).
#' 
#' If you are using exon microarray data, this step should be skipped, since exon 
#' microarray intensities are already on a continuous scale.
#' 
#' Previous step: Calculate normalization factors with \code{\link{runNormalization}}.
#' Next step: Initialize \code{\linkS4class{RegspliceResults}} object with the
#' constructor function \code{RegspliceResults()}.
#' 
#' 
#' @param rs_data \code{\linkS4class{RegspliceData}} object, which has been filtered with
#'   \code{\link{filterZeros}} and \code{\link{filterLowCounts}}, and (optionally) 
#'   normalization factors added with \code{\link{runNormalization}}.
#' 
#' 
#' @return Returns a \code{\linkS4class{RegspliceData}} object. Transformed counts are 
#'   stored in the \code{counts} matrix, and weights are stored in a new \code{weights} 
#'   data matrix. The data matrices can be accessed with the accessor functions 
#'   \code{\link{countsData}} and \code{\link{weightsData}}.
#' 
#' @seealso \code{\link{runNormalization}} \code{\link{fitRegMultiple}} 
#'   \code{\link{fitNullMultiple}} \code{\link{fitFullMultiple}}
#' 
#' @importFrom limma voom
#' @importFrom stats model.matrix
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
#' rs_data <- filterZeros(rs_data)
#' rs_data <- filterLowCounts(rs_data)
#' rs_data <- runNormalization(rs_data)
#' rs_data <- runVoom(rs_data)
#' 
runVoom <- function(rs_data) {
  
  norm_factors <- colData(rs_data)$norm_factors
  condition <- colData(rs_data)$condition
  counts <- countsData(rs_data)
  
  # design matrix
  design <- stats::model.matrix(~ condition)
  
  # normalized library sizes
  if (!is.null(norm_factors)) {
    lib_sizes <- norm_factors * colSums(counts)
  } else {
    lib_sizes <- NULL
  }
  
  colData(rs_data)$lib_sizes <- lib_sizes
  
  # run voom
  out_voom <- limma::voom(counts = counts, design = design, lib.size = lib_sizes)
  
  assays(rs_data)$counts <- out_voom$E
  assays(rs_data)$weights <- out_voom$weights
  
  rs_data
}



