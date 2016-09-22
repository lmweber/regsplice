library(regsplice)
context("Input data preparation")

test_that("error if number of rows in count table does not match total number of exons", {
  
  set.seed(123)
  counts <- matrix(sample(100:200, 7 * 4, replace = TRUE), nrow = 7, ncol = 4)
  n_exons <- c(3, 1, 2)
  gene_IDs <- paste0("gene", 1:length(n_exons))
  condition <- rep(c(0, 1), each = 2)
  
  expect_error(RegspliceData(counts, gene_IDs, n_exons, condition), "total number of exons")
})




test_that("genes containing only a single exon are removed", {
  
  set.seed(123)
  n_exons <- c(11, 1, 3, 1, 7)
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  gene_IDs <- paste0("gene", 1:length(n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  Y <- RegspliceData(counts, gene_IDs, n_exons, condition)
  Y <- filter_zeros(Y)
  Y <- filter_low_counts(Y)
  
  n_genes <- length(table(rowData(Y)$gene_IDs))
  n_exons_total <- nrow(rowData(Y))
  
  expect_equal(n_genes, 3)
  expect_equal(n_exons_total, 21)
})




test_that("genes with all zero counts are removed", {
  
  counts <- rbind(matrix(1, nrow = 5, ncol = 4), 
                  matrix(0, nrow = 3, ncol = 4), 
                  matrix(1, nrow = 2, ncol = 4))
  gene_IDs <- paste0("gene", 1:4)
  n_exons <- c(3, 2, 3, 2)
  condition <- rep(c(0, 1), each = 2)
  
  Y <- RegspliceData(counts, gene_IDs, n_exons, condition)
  Y <- filter_zeros(Y)
  
  n_genes <- length(table(rowData(Y)$gene_IDs))
  n_exons_total <- nrow(rowData(Y))
  
  expect_equal(n_genes, 3)
  expect_equal(n_exons_total, 7)
})



