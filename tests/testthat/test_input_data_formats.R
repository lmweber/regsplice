library(regsplice)
context("Input data formats")

test_that("counts matrix can be input as matrix or data frame", {
  
  counts <- matrix(1, nrow = 5, ncol = 6)
  counts_df <- as.data.frame(counts)
  gene_IDs <- c("gene1", "gene2")
  n_exons <- c(3, 2)
  condition <- rep(c(0, 1), each = 3)
  
  rs_data_matrix <- RegspliceData(counts,    gene_IDs, n_exons, condition)
  rs_data_df     <- RegspliceData(counts_df, gene_IDs, n_exons, condition)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = S4Vectors::SimpleList(counts = counts), 
    rowData = S4Vectors::DataFrame(gene_IDs = rep(gene_IDs, n_exons)), 
    colData = S4Vectors::DataFrame(condition))
  
  rs_se <- RegspliceData(se)
  
  expect_equivalent(rs_data_matrix, rs_data_df)
  expect_equivalent(rs_data_matrix, rs_se)
  expect_equivalent(rs_data_df, rs_se)
})



