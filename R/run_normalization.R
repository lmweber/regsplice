#' Calculate normalization factors.
#' 
#' Calculate normalization factors to scale library sizes, using the TMM (trimmed mean of
#' M-values) method implemented in \code{edgeR}.
#' 
#' Normalization factors are used to scale the raw library sizes (total read counts per 
#' sample). We use the TMM (trimmed mean of M-values) normalization method (Robinson and
#' Oshlack, 2010), as implemented in the \code{edgeR} package.
#' 
#' For more details, see the documentation for \code{\link[edgeR]{calcNormFactors}} in
#' the \code{edgeR} package.
#' 
#' This step should be performed after filtering. The normalization factors are then
#' passed on to \code{limma-voom} in the next step (\code{\link{run_voom}}).
#' 
#' Normalization should be skipped when using exon microarray data. (When using the 
#' \code{\link{regsplice}} wrapper function, normalization can be disabled with the
#' argument \code{normalize = FALSE}).
#' 
#' 
#' @param Y RNA-seq read counts for multiple genes (list of data frames or matrices). 
#'   Names contain gene names. Created using \code{\link{prepare_data}} and
#'   \code{\link{filter_exons}}.
#' @param norm_method Normalization method to use. Options are \code{"TMM"}, 
#'   \code{"RLE"}, \code{"upperquartile"}, and \code{"none"}. See documentation for 
#'   \code{\link[edgeR]{calcNormFactors}} in \code{edgeR} package for details. Default is
#'   \code{"TMM"}.
#' 
#' @return Returns a numeric vector of normalization factors (one value per sample).
#' 
#' @seealso \code{\link{prepare_data}} \code{\link{filter_exons}}
#'   \code{\link{run_voom}}
#' 
#' @importFrom edgeR calcNormFactors
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
run_normalization <- function(Y, norm_method = "TMM") {
  
  norm_method <- match.arg(norm_method, c("TMM", "RLE", "upperquartile", "none"))
  
  # collapse data
  counts <- do.call(rbind, Y)
  rownames(counts) <- NULL
  
  # run normalization
  edgeR::calcNormFactors(counts, method = norm_method)
}


