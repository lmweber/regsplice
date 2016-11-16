library(regsplice)
context("when_null_selected options")

test_that('when_null_selected option "ones" works correctly', {
  
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  set.seed(123)
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  gene_IDs <- paste0("gene", 1:length(n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  rs_data <- RegspliceData(counts, gene_IDs, n_exons, condition)
  
  rs_results <- regsplice(rs_data, when_null_selected = "ones")
  
  
  n <- 6  # expect 6 remaining genes
  
  expect_length(rs_results@gene_IDs, n)
  
  expect_length(rs_results@p_vals, n)
  expect_length(rs_results@p_adj, n)
  
  expect_true(sum(rs_results@p_vals == 1) > 0)
  expect_true(sum(rs_results@p_adj == 1) > 0)
})




test_that('when_null_selected option "full" works correctly', {
  
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  set.seed(123)
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  gene_IDs <- paste0("gene", 1:length(n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  rs_data <- RegspliceData(counts, gene_IDs, n_exons, condition)
  
  rs_results <- regsplice(rs_data, when_null_selected = "full")
  
  
  n <- 6  # expect 6 remaining genes
  
  expect_length(rs_results@gene_IDs, n)
  
  expect_length(rs_results@p_vals, 6)
  expect_length(rs_results@p_adj, 6)
  
  expect_true(sum(is.na(rs_results@p_vals)) == 0)
  expect_true(sum(is.na(rs_results@p_adj)) == 0)
})




test_that('when_null_selected option "NA" works correctly', {
  
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  set.seed(123)
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  gene_IDs <- paste0("gene", 1:length(n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  rs_data <- RegspliceData(counts, gene_IDs, n_exons, condition)
  
  rs_results <- regsplice(rs_data, when_null_selected = "NA")
  
  
  n <- 6  # expect 6 remaining genes
  
  expect_length(rs_results@gene_IDs, n)
  
  expect_length(rs_results@p_vals, 6)
  expect_length(rs_results@p_adj, 6)
  
  expect_true(sum(is.na(rs_results@p_vals)) > 0)
  expect_true(sum(is.na(rs_results@p_adj)) > 0)
})



