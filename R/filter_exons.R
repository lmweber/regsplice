#' Filter low-count exons.
#' 
#' Filter low-count exons from RNA-seq read count data.
#' 
#' Filters low-count exons from an RNA-seq read count data set prepared with 
#' \code{\link{prepare_data}}, i.e. in the format required by other functions in the
#' \code{regsplice} pipeline.
#' 
#' The arguments \code{n1} and \code{n2} control the amount of filtering. Exons that meet
#' the filtering conditions specified by \code{n1} and \code{n2} are kept. Default values
#' are provided; however these may not be optimal for some experimental designs.
#' 
#' Any single-exon genes that remain after exons have been removed during the filtering 
#' step are also removed (since differential splicing requires multiple exons). The
#' output is in the same format as from \code{\link{prepare_data}}.
#' 
#' @param Y RNA-seq read counts for multiple genes (list of data frames or matrices).
#'   Names contain gene names. Created using \code{\link{prepare_data}}.
#' @param n1 Filtering parameter: minimum number of reads per exon, summed across all 
#'   biological samples. Default is 6.
#' @param n2 Filtering parameter: minimum number of reads for a single biological sample
#'   per exon, i.e. at least one sample must have this number of reads. Default is 3.
#' 
#' @return Returns a list of data frames, where each data frame in the list contains the 
#'   RNA-seq read counts for one gene. Gene names are stored as names of the list items.
#'   Low-count exons and any remaining single-exon genes have been removed.
#' 
#' @family prepare_data filter_exons voom_weights
#' 
#' @export
#' 
#' @examples 
#' file_counts <- system.file("extdata/vignette_counts.txt", package = "regsplice")
#' data <- read.table(file_counts, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#' counts <- data[, 2:7]
#' gene <- sapply(strsplit(data$exon, ":"), function(s) s[[1]])
#' 
#' Y <- prepare_data(counts, gene)
#' Y <- filter_exons(Y)
#' 
filter_exons <- function(Y, n1 = 6, n2 = 3) {
  
  if (!(is.list(Y) & !is.data.frame(Y))) {
    stop("data Y for multiple genes must be a list of data frames or matrices")
  }
  
  # keep vector of gene names
  gene <- names(Y)
  n_exons <- sapply(Y, nrow)
  gene <- rep(gene, times = n_exons)
  
  # collapse list of data frames
  counts <- do.call(rbind, Y)
  
  filt_exons_n1 <- rowSums(counts) >= n1
  filt_exons_n2 <- apply(counts, MARGIN = 1, FUN = function(r) max(r) >= n2)
  ix_keep <- filt_exons_n1 & filt_exons_n2
  
  counts <- counts[ix_keep, , drop = FALSE]
  row.names(counts) <- NULL
  gene <- gene[ix_keep]
  
  if (nrow(counts) != length(gene)) {
    stop("Number of rows in 'counts' does not match length of 'gene' after filtering.")
  }
  
  n_filtered <- sum(ix_keep == FALSE)
  message(paste("removed", n_filtered, "low-count exon(s)"))
  
  # split count table again (as in prepare_data)
  Y <- split_genes(counts, gene)
  
  # remove any new single-exon genes that remain after filtering
  Y <- filter_genes_single_exon(Y)
}


