library(regsplice)
context("Input data preparation and filtering")

test_that(paste0("internal function 'ix_exons_zero_counts' gives indices of exons (rows) ", 
                 "with zero counts"), {
  counts <- rbind(matrix(1, nrow = 2, ncol = 4), 
                  matrix(0, nrow = 1, ncol = 4), 
                  matrix(1, nrow = 2, ncol = 4), 
                  matrix(0, nrow = 3, ncol = 4), 
                  matrix(1, nrow = 2, ncol = 4))
  
  n_zeros <- sum(ix_exons_zero_counts(counts = counts))
  
  expect_equivalent(n_zeros, 4)
})



test_that(paste0("internal function 'split_genes' returns error if length of genes argument ", 
                 "does not match number of rows in count table"), {
  counts <- matrix(sample(100:200, 7 * 4, replace = TRUE), nrow = 7, ncol = 4)
  gene <- paste0("gene", rep(1:3, times = c(3, 1, 2)))
  
  expect_error(split_genes(counts, gene), "Length")
})



test_that(paste0("internal function 'filter_genes_single_exon' removes genes containing ", 
                 "only a single exon"), {
  n_exons <- c(11, 1, 3, 1, 7)
  n_genes <- length(n_exons)
  gene <- paste0("gene", rep(1:n_genes, times = n_exons))
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  
  Y <- split_genes(counts = counts, gene = gene)
  Y_filtered <- filter_genes_single_exon(Y = Y)
  
  expect_length(Y, 5)
  expect_length(Y_filtered, 3)
})



test_that("exported function 'prepare_data' combines data preparation steps correctly", {
  counts <- rbind(matrix(1, nrow = 2, ncol = 4), 
                  matrix(0, nrow = 1, ncol = 4), 
                  matrix(1, nrow = 1, ncol = 4), 
                  matrix(0, nrow = 5, ncol = 4), 
                  matrix(1, nrow = 4, ncol = 4))
  gene <- paste0("gene", rep(1:7, times = c(2, 1, 1, 3, 2, 2, 2)))
  
  ix_zeros <- ix_exons_zero_counts(counts = counts)
  counts <- counts[!ix_zeros, ]
  gene <- gene[!ix_zeros]
  
  Y <- split_genes(counts = counts, gene = gene)
  Y_filtered <- filter_genes_single_exon(Y = Y)
  
  zeros_none <- function(d) all(d != 0)
  
  expect_length(Y, 4)
  expect_length(Y_filtered, 3)
  
  expect_true(all(sapply(Y, zeros_none)))
  expect_true(all(sapply(Y_filtered, zeros_none)))
})



test_that("genes with all zero counts are removed", {
  counts <- rbind(matrix(1, nrow = 5, ncol = 4), 
                  matrix(0, nrow = 3, ncol = 4), 
                  matrix(1, nrow = 2, ncol = 4))
  gene <- paste0("gene", rep(1:4, times = c(3, 2, 3, 2)))
  
  Y <- prepare_data(counts = counts, gene = gene)
  
  expect_length(Y, 3)
})


