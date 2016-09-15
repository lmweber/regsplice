# Internal functions for data preparation. These functions should not be called directly 
# by users. Use the exported function "prepare_data" instead.


# Given an RNA-seq read count table (matrix or data frame), return a logical vector
# indicating exons (rows) with zero counts in all biological samples (columns).
# 
ix_exons_zero_counts <- function(counts) {
  
  ix_zeros <- apply(counts, MARGIN = 1, function(d) all(d == 0))
  
  # number of exons with zero counts
  n_zeros <- sum(ix_zeros)
  message(paste("removed", n_zeros, "exon(s) with zero counts"))
  
  ix_zeros
}



# Function to split a count table (matrix or data frame) into a list of sub-tables (data 
# frames), one for each gene. The "gene" argument is a vector of gene IDs, with one entry
# for every row (exon) in the count table; i.e. its length should equal the number of
# rows in the count table.
# 
split_genes <- function(counts, gene) {
  
  if(length(gene) != nrow(counts)) {
    stop("length of vector of gene IDs is not equal to number of rows in count table")
  }
  
  # use data frame instead of matrix for counts so each sub-matrix keeps its shape
  if (is.matrix(counts)) counts <- as.data.frame(counts)
  
  if (is.null(colnames(counts))) colnames(counts) <- paste0("sample", 1:ncol(counts))
  
  # set factor levels to keep genes in original order
  if (!is.character(gene)) gene <- as.character(gene)
  gene <- factor(gene, levels = unique(gene))
  
  Y <- split(counts, gene)
}



# Filter genes containing only a single exon, from a data set prepared with
# "split_genes". (Differential splicing requires multi-exon genes.)
# 
filter_genes_single_exon <- function(Y) {
  
  single_exons <- sapply(Y, function(d) nrow(d) == 1)
  Y <- Y[!single_exons]
  
  # number of single-exon genes (note: exons with zero counts should have already been
  # removed with 'ix_exons_zero_counts')
  n_single_exons <- sum(single_exons)
  message(paste("removed", n_single_exons, "remaining single-exon gene(s)"))
  
  Y
}


