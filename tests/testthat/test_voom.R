library(regsplice)
context("voom transformation and weights")

test_that("voom transformation and weights are calculated correctly", {
  
  counts <- cbind(sample1 = c(1, 4, 0, 5), 
                  sample2 = c(0, 1, 0, 5), 
                  sample3 = c(0, 9, 6, 30), 
                  sample4 = c(3, 12, 4, 23), 
                  sample5 = c(1, 4, 0, 10), 
                  sample6 = c(1, 3, 1, 22))
  gene_IDs <- "gene1"
  n_exons <- 4
  condition <- rep(c(0, 1), each = 3)
  
  Y <- RegspliceData(counts, gene_IDs, n_exons, condition)
  Y <- filter_zeros(Y)
  Y <- filter_low_counts(Y)
  Y <- run_normalization(Y)
  Y <- run_voom(Y)
  
  lib_sizes <- colData(Y)$lib_sizes
  
  
  expect_length(lib_sizes, 6)
  
  expect_is(weightsData(Y), "matrix")
  expect_is(weightsData(Y)[2, 3], "numeric")
})



