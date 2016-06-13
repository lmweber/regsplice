#' Prepare data.
#' 
#' Prepare data into format required by \code{regsplice} functions.
#' 
#' Prepares data into the format required by other functions in the \code{regsplice}
#' pipeline.
#' 
#' Inputs are a matrix or data frame of RNA-seq read counts or exon microarray
#' intensities (\code{counts}) and a vector of gene IDs (\code{gene}). The vector
#' \code{gene} must have length equal to the number of rows in \code{counts}; i.e. one
#' entry for each exon, with repeated entries for multiple exons within the same gene.
#' The repeated entries are used to determine gene length.
#' 
#' Data preparation consists of the following steps:
#' \itemize{
#' \item Remove exons (rows) with zero counts in all biological samples (columns).
#' \item Split count or intensity table into a list of sub-tables (data frames), one for
#' each gene.
#' \item Remove genes containing only a single exon (since differential splicing requires
#' multiple exons). For this to work correctly, this step must occur after zero-count 
#' exons have already been removed.}
#' 
#' The function \code{\link{filter_exons}} can then be used to filter low-count exons.
#' 
#' @param counts RNA-seq read counts or exon microarray intensities (matrix or data
#'   frame). Each row is an exon, and each column is a biological sample.
#' @param gene Vector of gene IDs (character vector). Length is equal to the number of 
#'   rows in \code{counts}; i.e. one entry for each exon, with repeated entries for
#'   multiple exons within the same gene. The repeated entries are used to determine gene
#'   length.
#' 
#' @return Returns a list of data frames, where each data frame in the list contains the 
#'   RNA-seq read counts or exon microarray intensities for one gene. Gene names are
#'   stored as names of the list items. Exons with zero counts and single-exon genes have
#'   been removed.
#' 
#' @seealso \code{\link{filter_exons}} \code{\link{voom_weights}}
#' 
#' @export
#' 
#' @examples
#' # ---------
#' # Example 1
#' # ---------
#' 
#' counts <- matrix(sample(100:200, 14 * 4, replace = TRUE), nrow = 14, ncol = 4)
#' counts[c(2, 3, 7), ] <- 0
#' gene <- paste0("gene", rep(1:5, times = c(3, 2, 3, 1, 5)))
#' 
#' # show raw data
#' data.frame(gene = gene, counts)
#' 
#' Y <- prepare_data(counts, gene)
#' 
#' # show prepared data and number of genes
#' Y
#' length(Y)
#' 
#' 
#' # --------------------
#' # Example 2 (Vignette)
#' # --------------------
#' 
#' file_counts <- system.file("extdata/vignette_counts.txt", package = "regsplice")
#' data <- read.table(file_counts, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#' counts <- data[, 2:7]
#' gene <- sapply(strsplit(data$exon, ":"), function(s) s[[1]])
#' 
#' Y <- prepare_data(counts, gene)
#' 
prepare_data <- function(counts, gene) {
  
  # remove exons (rows in count table) with zero counts in all samples
  ix_zeros <- ix_exons_zero_counts(counts)
  counts <- counts[!ix_zeros, ]
  gene <- gene[!ix_zeros]
  
  # split count table into list of sub-tables, one for each gene
  Y <- split_genes(counts, gene)
  
  # remove single-exon genes (note the sequence of steps here, i.e. we remove any
  # single-exon genes that remain after removing exons with zero counts)
  Y <- filter_genes_single_exon(Y)
}


