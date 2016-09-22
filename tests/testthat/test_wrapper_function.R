library(regsplice)
context("regsplice wrapper function")

test_that("regsplice wrapper function gives correct results", {
  
  set.seed(123)
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  gene_IDs <- paste0("gene", 1:length(n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  res <- regsplice(counts, gene_IDs, n_exons, condition)
  
  
  n <- 6  # expect 6 remaining genes
  
  expect_length(res@gene_IDs, n)
  
  expect_length(res@p_val, n)
  expect_length(res@p_adj, n)
  expect_length(res@LR_stat, n)
  expect_length(res@df_test, n)
  
  expect_true(all(!is.na(res@p_val)))
  expect_true(all(!is.na(res@p_adj)))
  
  expect_true(all(res@p_val >= 0))
  expect_true(all(res@p_val <= 1))
  expect_true(all(res@p_adj >= 0))
  expect_true(all(res@p_adj <= 1))
  
  expect_true(all(res@LR_stat >= 0, na.rm = TRUE))
  expect_true(all(res@df_test >= 1, na.rm = TRUE))
})



