library(regsplice)
context("Normalization factors")

test_that("normalization factors work correctly", {
  
  counts <- rbind(c(101, 198, 53, 102, 100, 101), 
                  c(90, 220, 60, 89, 90, 103), 
                  c(95, 205, 66, 105, 85, 99), 
                  c(300, 622, 145, 311, 290, 320), 
                  c(290, 602, 140, 290, 303, 298), 
                  c(290, 410, 90, 190, 211, 209))
  gene_IDs <- paste0("gene", 1:2)
  n_exons <- c(3, 3)
  condition <- rep(c(0, 1), each = 3)
  
  Y <- RegspliceData(counts, gene_IDs, n_exons, condition)
  Y <- filter_zeros(Y)
  Y <- filter_low_counts(Y)
  Y <- run_normalization(Y)
  
  norm_factors <- colData(Y)$norm_factors
  
  expect_length(norm_factors, 6)
  expect_is(norm_factors, "numeric")
})



