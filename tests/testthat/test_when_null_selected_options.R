library(regsplice)
context("when_null_selected options")

test_that('when_null_selected option "ones" works correctly', {
  
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  set.seed(123)
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  gene_IDs <- paste0("gene", 1:length(n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  res <- regsplice(counts, gene_IDs, n_exons, condition, when_null_selected = "ones")
  
  
  n <- 6  # expect 6 remaining genes
  
  expect_length(res@gene_IDs, n)
  
  expect_length(res@p_val, n)
  expect_length(res@p_adj, n)
  
  expect_true(sum(res@p_val == 1) > 0)
  expect_true(sum(res@p_adj == 1) > 0)
})




test_that('when_null_selected option "full" works correctly', {
  
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  set.seed(123)
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  gene_IDs <- paste0("gene", 1:length(n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  res <- regsplice(counts, gene_IDs, n_exons, condition, when_null_selected = "full")
  
  
  n <- 6  # expect 6 remaining genes
  
  expect_length(res@gene_IDs, n)
  
  expect_length(res@p_val, 6)
  expect_length(res@p_adj, 6)
  
  expect_true(sum(is.na(res@p_val)) == 0)
  expect_true(sum(is.na(res@p_adj)) == 0)
})




test_that('when_null_selected option "NA" works correctly', {
  
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  set.seed(123)
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  gene_IDs <- paste0("gene", 1:length(n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  res <- regsplice(counts, gene_IDs, n_exons, condition, when_null_selected = "NA")
  
  
  n <- 6  # expect 6 remaining genes
  
  expect_length(res@gene_IDs, n)
  
  expect_length(res@p_val, 6)
  expect_length(res@p_adj, 6)
  
  expect_true(sum(is.na(res@p_val)) > 0)
  expect_true(sum(is.na(res@p_adj)) > 0)
})



