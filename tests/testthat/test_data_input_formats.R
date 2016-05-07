library(regsplice)
context("Input data loads correctly")

test_that("counts matrix can be input as matrix or data frame", {
  counts <- matrix(1, nrow = 5, ncol = 4)
  counts_df <- as.data.frame(counts)
  gene <- rep(c("gene1", "gene2"), times = c(3, 2))
  
  Y_matrix <- split_genes(counts = counts,    gene = gene)
  Y_df     <- split_genes(counts = counts_df, gene = gene)
  
  expect_equal(Y_matrix, Y_df)
})


test_that("gene identifiers can be input in several ways", {
  counts <- matrix(1, nrow = 5, ncol = 4)
  
  gene_chr <- rep(c("gene1", "gene2"), times = c(3, 2))
  gene_num <- rep(c(0, 1), times = c(3, 2))
  gene_fac <- as.factor(rep(c(0, 1), times = c(3, 2)))
  
  Y_chr <- split_genes(counts = counts, gene = gene_chr)
  Y_num <- split_genes(counts = counts, gene = gene_num)
  Y_fac <- split_genes(counts = counts, gene = gene_fac)
  
  expect_equal(unname(Y_chr), unname(Y_num))
  expect_equal(unname(Y_chr), unname(Y_fac))
})

