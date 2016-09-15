library(regsplice)
context("when_null_selected options")

test_that('when_null_selected option "ones" works correctly', {
  
  # generate data
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  gene <- paste0("gene", rep(1:length(n_exons), times = n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  set.seed(123)
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  
  res <- regsplice(counts = counts, gene = gene, condition = condition, 
                   when_null_selected = "ones")
  
  n_genes <- 6
  
  expect_length(res$p_val, 6)
  expect_length(res$p_adj, 6)
  
  expect_true(sum(res$p_val == 1) > 0)
  expect_true(sum(res$p_adj == 1) > 0)
})



test_that('when_null_selected option "GLM" works correctly', {
  
  # generate data
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  gene <- paste0("gene", rep(1:length(n_exons), times = n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  set.seed(123)
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  
  res <- regsplice(counts = counts, gene = gene, condition = condition, 
                   when_null_selected = "GLM")
  
  n_genes <- 6
  
  expect_length(res$p_val, 6)
  expect_length(res$p_adj, 6)
  
  expect_true(sum(is.na(res$p_val)) == 0)
  expect_true(sum(is.na(res$p_adj)) == 0)
})



test_that('when_null_selected option "NA" works correctly', {
  
  # generate data
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  gene <- paste0("gene", rep(1:length(n_exons), times = n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  set.seed(123)
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  
  res <- regsplice(counts = counts, gene = gene, condition = condition, 
                   when_null_selected = "NA")
  
  n_genes <- 6
  
  expect_length(res$p_val, 6)
  expect_length(res$p_adj, 6)
  
  expect_true(sum(is.na(res$p_val)) > 0)
  expect_true(sum(is.na(res$p_adj)) > 0)
})


