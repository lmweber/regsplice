library(regsplice)
context("Filtering functions")

test_that("genes with zero counts are removed", {
  counts <- rbind(matrix(1, nrow = 5, ncol = 4), 
                  matrix(0, nrow = 3, ncol = 4), 
                  matrix(1, nrow = 2, ncol = 4))
  gene <- paste0("gene", rep(1:4, times = c(3, 2, 3, 2)))
  
  Y <- split_genes(counts = counts, gene = gene)
  Y_filtered <- filter_zeros(Y = Y)
  
  expect_length(Y, 4)
  expect_length(Y_filtered, 3)
})


test_that("genes containing only a single exon are removed", {
  n_exons <- c(11, 1, 3, 1, 7)
  n_genes <- length(n_exons)
  gene <- paste0("gene", rep(1:n_genes, times = n_exons))
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  
  Y <- split_genes(counts = counts, gene = gene)
  Y_filtered <- filter_single_exons(Y = Y)
  
  expect_length(Y, 5)
  expect_length(Y_filtered, 3)
})

