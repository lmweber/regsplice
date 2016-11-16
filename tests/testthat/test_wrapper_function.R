library(regsplice)
context("regsplice wrapper function")

test_that("regsplice wrapper function gives correct results", {
  
  # generate data
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  gene <- paste0("gene", rep(1:length(n_exons), times = n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  
  res <- regsplice(counts, gene, condition)
  
  
  n_genes <- 6  # non-single-exon genes only
  
  expect_length(res$p_vals, n_genes)
  expect_length(res$p_adj, n_genes)
  expect_length(res$LR_stats, n_genes)
  expect_length(res$df_tests, n_genes)
  
  expect_true(all(!is.na(res$p_vals)))
  expect_true(all(!is.na(res$p_adj)))
  
  expect_true(all(res$p_vals >= 0))
  expect_true(all(res$p_vals <= 1))
  expect_true(all(res$p_adj >= 0))
  expect_true(all(res$p_adj <= 1))
  
  expect_true(all(res$LR_stats >= 0, na.rm = TRUE))
  expect_true(all(res$df_tests >= 1, na.rm = TRUE))
})

