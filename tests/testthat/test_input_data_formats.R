library(regsplice)
context("Input data formats")

test_that("counts matrix can be input as matrix or data frame", {
  
  counts <- matrix(1, nrow = 5, ncol = 6)
  counts_df <- as.data.frame(counts)
  gene_IDs <- c("gene1", "gene2")
  n_exons <- c(3, 2)
  condition <- rep(c(0, 1), each = 3)
  
  Y_matrix <- RegspliceData(counts,    gene_IDs, n_exons, condition)
  Y_df     <- RegspliceData(counts_df, gene_IDs, n_exons, condition)
  
  expect_equivalent(Y_matrix, Y_df)
})



