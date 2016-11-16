#' @include class_RegspliceData.R class_RegspliceResults.R
NULL




#' Filter zero-count exons.
#' 
#' Filter exons with zero RNA-seq read counts in all biological samples.
#' 
#' Removes exon bins with zero RNA-seq read counts in all biological samples. Any 
#' remaining single-exon genes (after filtering) are also removed (since differential
#' splicing requires multiple exon bins).
#' 
#' Input data is assumed to be in the form of a \code{RegspliceData} object. See 
#' \code{\link{RegspliceData}} for details.
#' 
#' After filtering zero-count exon bins, any remaining genes containing only a single
#' exon bin are also removed (since differential splicing requires multiple exon bins).
#' 
#' Filtering should be skipped when using exon microarray data. (When using the 
#' \code{regsplice} wrapper function, filtering can be disabled with the argument 
#' \code{filter = FALSE}).
#' 
#' Previous step: Create \code{RegspliceData} object with \code{\link{RegspliceData}}
#' constructor function.
#' Next step: Filter low-count exon bins with \code{\link{filter_low_counts}}.
#' 
#' 
#' @param data \code{\linkS4class{RegspliceData}} object.
#' 
#' 
#' @return Returns a \code{\linkS4class{RegspliceData}} object.
#' 
#' @seealso \code{\linkS4class{RegspliceData}} \code{\link{filter_low_counts}}
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
#' Y <- RegspliceData(counts, gene_IDs, n_exons, condition)
#' 
#' Y <- filter_zeros(Y)
#' 
filter_zeros <- function(data) {
  
  if (!("RegspliceData" %in% is(data))) stop("'data' must be a 'RegspliceData' object")
  
  # remove exon bins (rows) with zero counts in all samples (columns)
  counts <- countsData(data)
  ix_zeros <- apply(counts, MARGIN = 1, function(d) all(d == 0))
  
  message(paste("removed", sum(ix_zeros), "exon(s) with zero counts"))
  
  data <- suppressMessages(data[!ix_zeros, ])
  
  # remove any remaining single-exon genes after filtering
  remove_single_exon_genes(data)
}



