library(regsplice)
context("Filter low-count exons")

test_that("low-count exon bins are filtered correctly", {
  
  counts <- rbind(c(1, 0, 0, 0, 0, 0), 
                  c(2, 3, 5, 1, 2, 1), 
                  c(0, 2, 1, 1, 2, 2), 
                  c(0, 4, 0, 0, 0, 0), 
                  c(0, 3, 4, 1, 1, 2), 
                  c(10, 5, 6, 8, 2, 8), 
                  c(13, 8, 9, 2, 4, 2))
  gene_IDs <- paste0("gene", 1:3)
  n_exons <- c(3, 2, 2)
  condition <- rep(c(0, 1), each = 3)
  
  Y <- RegspliceData(counts, gene_IDs, n_exons, condition)
  
  Y <- filter_zeros(Y)
  Y <- filter_low_counts(Y)
  
  n_genes <- length(names(table(rowData(Y)$gene_IDs)))
  n_exons <- nrow(rowData(Y))
  
  expect_equivalent(n_genes, 1)
  expect_equivalent(n_exons, 2)
})



