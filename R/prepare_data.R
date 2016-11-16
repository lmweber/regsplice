#' Split data into list of genes.
#' 
#' Function to prepare input data by splitting an RNA-seq count table into a list of 
#' sub-tables, one for each gene. This is the format required by the \emph{regsplice} 
#' pipeline.
#' 
#' Splits a matrix (or data frame) of RNA-seq counts into a list of data frames, where 
#' each data frame in the list contains the RNA-seq counts for one gene.
#' 
#' Gene IDs are provided in the \code{gene} argument, which should contain one entry for
#' every exon in every gene; i.e. the length of \code{gene} should be equal to the number
#' of rows in \code{counts}, with repeated entries for for multiple exons in each gene.
#'
#' @param counts RNA-seq counts (matrix or data frame). Each row is an exon, and each
#'   column is a biological sample.
#' @param gene Vector of gene IDs (character vector). Length is equal to the number of
#'   rows in \code{counts}.
#'
#' @return Returns a list of data frames, where each data frame contains the RNA-seq 
#'   counts for one gene.
#' 
#' @family split_genes filter_zeros
#' 
#' @export
#'
#' @examples
#' counts <- matrix(sample(100:200, 7 * 4, replace = TRUE), nrow = 7, ncol = 4)
#' gene <- paste0("gene", rep(1:3, times = c(3, 2, 2)))
#' Y <- split_genes(counts, gene)
#' 
#' file_counts <- system.file("extdata/counts.txt", package = "regsplice")
#' data <- read.table(file_counts, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#' counts <- data[, 2:7]
#' gene <- sapply(strsplit(data$exon, ":"), function(s) s[[1]])
#' Y <- split_genes(counts, gene)
#' 
split_genes <- function(counts, gene) {
  
  # use data frame instead of matrix for counts so each sub-matrix keeps its shape
  if (is.matrix(counts)) counts <- as.data.frame(counts)
  if (is.null(colnames(counts))) colnames(counts) <- paste0("sample", 1:ncol(counts))
  
  # set factor levels to keep genes in original order
  if (!is.character(gene)) gene <- as.character(gene)
  gene <- factor(gene, levels = unique(gene))
  
  Y <- split(counts, gene)
}



#' Filter genes with zero counts.
#' 
#' Filter genes with zero counts from a data set prepared with \code{\link{split_genes}}.
#' 
#' Input data should be in the format prepared by the \code{\link{split_genes}} function,
#' i.e. a list of data frames, where each data frame contains RNA-seq counts for one
#' gene.
#'
#' @param Y RNA-seq count table in the format prepared by \code{\link{split_genes}}; i.e.
#'   a list of data frames, one for each gene.
#'
#' @return Y Returns a list of data frames, where each data frame contains the RNA-seq 
#'   counts for one gene, and genes with zero total counts have been removed.
#' 
#' @family split_genes filter_zeros filter_single_exons
#' 
#' @export
#'
#' @examples
#' counts <- rbind(matrix(1, nrow = 3, ncol = 4), 
#'                 matrix(1, nrow = 2, ncol = 4), 
#'                 matrix(0, nrow = 4, ncol = 4), 
#'                 matrix(1, nrow = 2, ncol = 4))
#' gene <- paste0("gene", rep(1:4, times = c(3, 2, 4, 2)))
#' Y <- split_genes(counts, gene)
#' names(Y)
#' length(Y)
#' Y <- filter_zeros(Y)
#' names(Y)
#' length(Y)
#' 
filter_zeros <- function(Y) {
  # genes with zero counts
  zeros <- sapply(Y, function(d) all(d == 0))
  Y <- Y[!zeros]
}



#' Filter genes containing a single exon.
#' 
#' Filter genes containing only a single exon, from a data set prepared with
#' \code{\link{split_genes}}.
#' 
#' Input data should be in the format prepared by the \code{\link{split_genes}} function,
#' i.e. a list of data frames, where each data frame contains RNA-seq counts for one
#' gene.
#'
#' @param Y RNA-seq count table in the format prepared by \code{\link{split_genes}}; i.e.
#'   a list of data frames, one for each gene.
#'
#' @return Y Returns a list of data frames, where each data frame contains the RNA-seq 
#'   counts for one gene, and genes containing only a single exon have been removed.
#' 
#' @family split_genes filter_zeros filter_single_exons
#' 
#' @export
#'
#' @examples
#' counts <- rbind(matrix(1, nrow = 3, ncol = 4), 
#'                 matrix(1, nrow = 2, ncol = 4), 
#'                 matrix(1, nrow = 1, ncol = 4), 
#'                 matrix(0, nrow = 5, ncol = 4), 
#'                 matrix(1, nrow = 2, ncol = 4))
#' gene <- paste0("gene", rep(1:5, times = c(3, 2, 1, 5, 2)))
#' Y <- split_genes(counts, gene)
#' names(Y)
#' length(Y)
#' Y <- filter_single_exons(Y)
#' names(Y)
#' length(Y)
#' 
filter_single_exons <- function(Y) {
  # genes containing a single exon
  single_exons <- sapply(Y, function(d) nrow(d) == 1)
  Y <- Y[!single_exons]
}


