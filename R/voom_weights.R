#' Calculate 'voom' weights, transformation, and normalization.
#' 
#' Calculate exon-level precision weights, log2-counts per million transformation, and 
#' scale normalization across samples, using \code{limma-voom}.
#' 
#' In many cases, power to detect differential exon usage can be further improved by 
#' using precision weights and normalizing across samples. This function uses
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
#' Note that \code{voom} assumes that exons (rows) with zero or low counts have already 
#' been removed, so this function should be used after data preparation and filtering 
#' with \code{\link{prepare_data}} and \code{\link{filter_exons}}.
#' 
#' If you are running \code{regsplice} with the wrapper function \code{\link{regsplice}},
#' this function will be called automatically, and will perform all three calculations 
#' (weights, continuous transformation, normalization). However, depending on your data, 
#' you may wish to use only the weights. This can be done by running the individual 
#' functions in the \code{regsplice} workflow; after running \code{voom_weights}, extract
#' the \code{weights} object from the returned list object, and pass it to the model 
#' fitting functions using the \code{weights} arguments. See the example workflow in the
#' vignette for more details.
#' 
#' 
#' @param Y RNA-seq read counts for multiple genes (list of data frames or matrices), 
#'   after preparation and filtering with \code{\link{prepare_data}} and
#'   \code{\link{filter_exons}}. Note that \code{voom} assumes filtered data (see above).
#' @param condition Experimental conditions for each sample (character or numeric vector,
#'   or factor).
#' 
#' @return Returns a list containing:
#' \itemize{
#' \item Y: Transformed and normalized data (log2-counts per million; normalized using 
#' scale normalization). List of data frames, where each data frame contains the data for
#' one gene.
#' \item weights: Exon-level precision weights. List of data frames, where each data
#' frame contains the weights for one gene.
#' }
#' 
#' @family prepare_data filter_exons voom_weights
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
voom_weights <- function(Y, condition) {
  
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
  
  list(Y = out_counts, weights = out_weights)
}


