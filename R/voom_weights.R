#' Calculate 'voom' weights, transformation, and normalization.
#' 
#' Calculate exon-level precision weights, as well as optional log2-counts per million 
#' continuous transformation and scale normalization across samples; using
#' \code{limma-voom}.
#' 
#' In many cases, power to detect differential exon usage can be further improved by 
#' using precision weights and/or normalizing across samples. This function uses
#' \code{limma-voom} to do the following:
#' 
#' \itemize{
#' \item calculate exon-level precision weights
#' \item transform count data to log2-counts per million, i.e. on a continuous scale
#' \item normalize across samples using scale normalization
#' }
#' 
#' For more details, see the documentation for \code{\link[limma]{voom}} in the 
#' \code{limma} package.
#' 
#' By default, only the weights are returned. This is because \code{voom} calculates 
#' normalization factors on the continuous transformed data, so it is not possible to 
#' return normalization factors without transforming the data. The continuous 
#' transformation may be risky for some data sets, since it is a major modification of 
#' the data; so by default we do not use the transformed/normalized data. If you wish to 
#' use the transformed/normalized data, set the argument \code{norm = TRUE}.
#' 
#' Note that \code{voom} assumes that exons (rows) with zero or low counts have already 
#' been removed, so this function should be used after data preparation and filtering 
#' with \code{\link{prepare_data}} and \code{\link{filter_exons}}.
#' 
#' Exon microarray intensities should be log2-transformed prior to model fitting; this 
#' can either be done externally (for example with \code{voom}), or by setting \code{norm
#' = TRUE}.
#' 
#' 
#' @param Y RNA-seq read counts or exon microarray intensities for multiple genes (list
#'   of data frames or matrices; names contain gene names), after preparation and
#'   filtering with \code{\link{prepare_data}} and \code{\link{filter_exons}}. Note that
#'   \code{voom} assumes filtered data (see above).
#' @param condition Experimental conditions for each sample (character or numeric vector,
#'   or factor).
#' @param norm Whether to return continuous transformed/normalized data (log2-counts per
#'   million transformation; scale normalization across samples). Default is FALSE.
#'   Should be set to TRUE if using exon microarray intensities that have not already
#'   been log2-transformed.
#' 
#' @return Returns a list containing:
#' \itemize{
#' \item Y: RNA-seq read counts or exon microarray intensities (same as input data), or
#' continuous transformed/normalized data (log2-counts per million transformation; scale
#' normalization across samples). List of data frames, where each data frame contains the
#' data for one gene. Gene names are stored as names of the list items.
#' \item weights: Exon-level precision weights. List of data frames, where each data
#' frame contains the weights for one gene.
#' }
#' 
#' @seealso \code{\link{prepare_data}} \code{\link{filter_exons}}
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
#' condition <- rep(c(0, 1), each = 3)
#' 
#' Y <- prepare_data(counts, gene)
#' Y <- filter_exons(Y)
#' out_voom <- voom_weights(Y, condition)
#' 
voom_weights <- function(Y, condition, norm = FALSE) {
  
  # get gene names and collapse data frames
  n_exons <- sapply(Y, nrow)
  gene <- rep(names(Y), times = n_exons)
  counts <- do.call(rbind, Y)
  rownames(counts) <- NULL
  
  # run voom
  out_voom <- limma::voom(counts = counts, 
                          design = stats::model.matrix(~ condition), 
                          normalize.method = "scale")
  
  # format for regsplice functions
  out_counts <- split_genes(counts = out_voom$E, gene = gene)
  out_weights <- split_genes(counts = out_voom$weights, gene = gene)
  
  if (norm) {
    return(list(Y = out_counts, weights = out_weights))
  } else {
    return(list(Y = Y, weights = out_weights))
  }
}


