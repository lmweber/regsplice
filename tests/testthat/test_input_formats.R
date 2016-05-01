library(regsplice)
context("Input data loads correctly")

test_that("counts can be input as matrix or data frame", {
  counts <- matrix(1, nrow = 5, ncol = 4)
  counts_df <- as.data.frame(counts)
  gene <- rep(c("gene1", "gene2"), times = c(3, 2))
  condition <- rep(c(0, 1), each = 2)
  
  prepare_matrix <- prepare_data(counts = counts,    gene = gene, condition = condition)
  prepare_df     <- prepare_data(counts = counts_df, gene = gene, condition = condition)
  
  expect_s4_class(prepare_matrix, "SummarizedExperiment0")
  expect_s4_class(prepare_df,     "SummarizedExperiment0")
})


test_that("gene identifiers can be input in several ways", {
  counts <- matrix(1, nrow = 5, ncol = 4)
  condition <- rep(c(0, 1), each = 2)
  
  gene_chr <- rep(c("gene1", "gene2"), times = c(3, 2))
  gene_num <- rep(c(0, 1), times = c(3, 2))
  gene_fac <- as.factor(rep(c(0, 1), times = c(3, 2)))
  
  prepare_chr <- prepare_data(counts = counts, gene = gene_chr, condition = condition)
  prepare_num <- prepare_data(counts = counts, gene = gene_num, condition = condition)
  prepare_fac <- prepare_data(counts = counts, gene = gene_fac, condition = condition)
  
  expect_s4_class(prepare_chr, "SummarizedExperiment0")
  expect_s4_class(prepare_num, "SummarizedExperiment0")
  expect_s4_class(prepare_fac, "SummarizedExperiment0")
})


test_that("condition vector can be input in several ways", {
  counts <- matrix(1, nrow = 5, ncol = 4)
  gene <- rep(c("gene1", "gene2"), times = c(3, 2))
  
  condition_chr <- rep(c("condition0", "condition1"), each = 2)
  condition_num <- rep(c(0, 1), each = 2)
  condition_fac <- as.factor(rep(c(0, 1), each = 2))
  
  prepare_chr <- prepare_data(counts = counts, gene = gene, condition = condition_chr)
  prepare_num <- prepare_data(counts = counts, gene = gene, condition = condition_num)
  prepare_fac <- prepare_data(counts = counts, gene = gene, condition = condition_fac)
  
  expect_s4_class(prepare_chr, "SummarizedExperiment0")
  expect_s4_class(prepare_num, "SummarizedExperiment0")
  expect_s4_class(prepare_fac, "SummarizedExperiment0")
})


