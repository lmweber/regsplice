#' @include class_RegspliceData.R class_RegspliceResults.R
NULL




#' Filter low-count exons.
#' 
#' Filter low-count exons from RNA-seq read count data.
#' 
#' Filters low-count exon bins from RNA-seq read count data. Any remaining single-exon 
#' genes (after filtering) are also removed (since differential splicing requires
#' multiple exon bins).
#' 
#' Input data is assumed to be in the form of a \code{RegspliceData} object. See 
#' \code{\link{RegspliceData}} for details.
#' 
#' The arguments \code{filter_min_per_exon} and \code{filter_min_per_sample} control the 
#' amount of filtering. Exon bins that meet the filtering conditions are kept. Default 
#' values for the arguments are provided; however, these should be adjusted depending on 
#' the total number of samples and the number of samples per condition.
#' 
#' After filtering low-count exon bins, any remaining genes containing only a single exon
#' bin are also removed (since differential splicing requires multiple exon bins).
#' 
#' Filtering should be skipped when using exon microarray data. (When using the 
#' \code{regsplice} wrapper function, filtering can be disabled with the argument 
#' \code{filter = FALSE}).
#' 
#' Previous step: Filter zero-count exon bins with \code{\link{filterZeros}}.
#' Next step: Calculate normalization factors with \code{\link{runNormalization}}.
#' 
#' 
#' @param rs_data \code{\linkS4class{RegspliceData}} object.
#' @param filter_min_per_exon Filtering parameter: minimum number of reads per exon bin, 
#'   summed across all biological samples. Default is 6.
#' @param filter_min_per_sample Filtering parameter: minimum number of reads per 
#'   biological sample; i.e. for each exon bin, at least one sample must have this number
#'   of reads. Default is 3.
#' 
#' 
#' @return Returns a \code{\linkS4class{RegspliceData}} object.
#' 
#' @seealso \code{\link{filterZeros}} \code{\link{runNormalization}}
#' 
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
#' rs_data <- filterZeros(rs_data)
#' rs_data <- filterLowCounts(rs_data)
#' 
filterLowCounts <- function(rs_data, filter_min_per_exon = 6, filter_min_per_sample = 3) {
  
  if (!("RegspliceData" %in% is(rs_data))) stop("'rs_data' must be a 'RegspliceData' object")
  
  counts <- countsData(rs_data)
  filt_per_exon <- rowSums(counts) >= filter_min_per_exon
  filt_per_sample <- apply(counts, MARGIN = 1, FUN = function(r) max(r) >= filter_min_per_sample)
  
  ix_keep <- filt_per_exon & filt_per_sample
  
  message(paste("removed", sum(!ix_keep), "low-count exon(s)"))
  
  rs_data <- suppressMessages(rs_data[ix_keep, ])
  
  # remove any remaining single-exon genes after filtering
  .removeSingleExonGenes(rs_data)
}



