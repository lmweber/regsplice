library(regsplice)
context("regsplice wrapper function")

test_that("'regsplice()' wrapper function gives correct results", {
  
  set.seed(123)
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  gene_IDs <- paste0("gene", 1:length(n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  rs_data <- RegspliceData(counts, gene_IDs, n_exons, condition)
  
  rs_results <- regsplice(rs_data)
  
  
  n <- 6  # expect 6 remaining genes
  
  expect_length(rs_results@gene_IDs, n)
  
  expect_length(rs_results@p_vals, n)
  expect_length(rs_results@p_adj, n)
  expect_length(rs_results@LR_stats, n)
  expect_length(rs_results@df_tests, n)
  
  expect_true(all(!is.na(rs_results@p_vals)))
  expect_true(all(!is.na(rs_results@p_adj)))
  
  expect_true(all(rs_results@p_vals >= 0))
  expect_true(all(rs_results@p_vals <= 1))
  expect_true(all(rs_results@p_adj >= 0))
  expect_true(all(rs_results@p_adj <= 1))
  
  expect_true(all(rs_results@LR_stats >= 0, na.rm = TRUE))
  expect_true(all(rs_results@df_tests >= 1, na.rm = TRUE))
})



