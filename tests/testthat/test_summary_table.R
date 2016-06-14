library(regsplice)
context("Summary table")

test_that("summary table function works correctly", {
  
  # generate data
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  gene <- paste0("gene", rep(1:length(n_exons), times = n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  
  res <- regsplice(counts = counts, gene = gene, condition = condition)
  
  table_default             <- summary_table(res)
  table_all_up_to_threshold <- summary_table(res, n = Inf)
  table_all                 <- summary_table(res, n = Inf, threshold = 1)
  
  table_pvals <- summary_table(res, rank_by = "p-value")
  table_none  <- summary_table(res, rank_by = "none")
  
  table_small_default <- summary_table(res, n = 3)
  table_small_pvals   <- summary_table(res, n = 3, rank_by = "p-value")
  table_small_none    <- summary_table(res, n = 3, rank_by = "none")
  
  
  expect_equal(ncol(table_default), 5)
  expect_equal(ncol(table_all_up_to_threshold), 5)
  expect_equal(nrow(table_all), 6)
  expect_equal(ncol(table_all), 5)
  
  expect_equal(ncol(table_pvals), 5)
  expect_equal(nrow(table_none), 6)
  expect_equal(ncol(table_none), 5)
  
  expect_equal(ncol(table_small_default), 5)
  expect_equal(ncol(table_small_pvals), 5)
  expect_equal(nrow(table_small_none), 3)
  expect_equal(ncol(table_small_none), 5)
})

