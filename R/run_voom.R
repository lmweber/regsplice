#' Calculate 'voom' transformation and weights.
#' 
#' Use \code{limma-voom} to transform counts and calculate exon-level weights.
#' 
#' Raw integer counts do not fulfill the statistical assumptions required for linear 
#' modeling. The \code{limma-voom} methodology transforms counts to log2-counts per 
#' million (logCPM), and calculates exon-level weights based on the observed 
#' mean-variance relationship. Linear modeling methods can then be used with the 
#' transformed data.
#' 
#' For more details, see the documentation for \code{\link[limma]{voom}} in the 
#' \code{limma} package.
#' 
#' Note that \code{voom} assumes that exons (rows) with zero or low counts have already 
#' been removed, so this step should be done after filtering with
#' \code{\link{filter_exons}}.
#' 
#' Normalization factors from \code{\link{run_normalization}} can be provided with the 
#' \code{norm_factors} argument. These are used by \code{voom} to calculate normalized
#' library sizes. If the are not provided, \code{voom} will use non-normalized library
#' sizes (columnwise total counts) instead.
#' 
#' If you are using exon microarray data, this step should be skipped, since exon
#' microarray intensities are already on a continuous scale.
#' 
#' 
#' @param Y RNA-seq read counts for multiple genes (list of data frames or matrices;
#'   names contain gene names), after preparation and filtering with
#'   \code{\link{prepare_data}} and \code{\link{filter_exons}}. Note that \code{voom}
#'   assumes filtered data (see above).
#' @param condition Experimental conditions for each sample (character or numeric vector,
#'   or factor).
#' @param norm_factors Normalization factors from \code{\link{run_normalization}} 
#'   (numeric vector), which are used to calculate normalized library sizes. If not
#'   provided, non-normalized library sizes are used instead. Default is NULL.
#' 
#' @return Returns a list containing:
#' \itemize{
#' \item Y: Transformed count data (log2-counts per million, logCPM). List of data
#' frames, where each data frame contains the data for one gene. Gene names are stored as
#' names of the list items.
#' \item weights: Exon-level weights. List of data frames, where each data frame contains
#' the weights for one gene.
#' }
#' 
#' @seealso \code{\link{prepare_data}} \code{\link{filter_exons}}
#'   \code{\link{run_normalization}}
#' 
#' @importFrom limma voom
#' @importFrom stats model.matrix
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
#' Y <- prepare_data(counts, gene)
#' Y <- filter_exons(Y)
#' norm_factors <- run_normalization(Y)
#' 
#' out_voom <- run_voom(Y, condition, norm_factors = norm_factors)
#' Y <- out_voom$Y
#' weights <- out_voom$weights
#' 
run_voom <- function(Y, condition, norm_factors = NULL) {
  
  # get gene names and collapse data frames
  n_exons <- sapply(Y, nrow)
  gene <- rep(names(Y), times = n_exons)
  counts <- do.call(rbind, Y)
  rownames(counts) <- NULL
  
  # design matrix
  design <- stats::model.matrix(~ condition)
  
  # normalized library sizes
  if (!is.null(norm_factors)) {
    lib_sizes <- norm_factors * colSums(counts)
  } else {
    lib_sizes <- NULL
  }
  
  # run voom
  out_voom <- limma::voom(counts = counts, design = design, lib.size = lib_sizes)
  
  # format for regsplice functions
  out_counts <- split_genes(counts = out_voom$E, gene = gene)
  out_weights <- split_genes(counts = out_voom$weights, gene = gene)
  
  list(Y = out_counts, weights = out_weights)
}


