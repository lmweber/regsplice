library(regsplice)
context("Vignette example")

test_that("results from vignette example are as expected", {
  
  file_counts <- system.file("extdata/vignette_counts.txt", package = "regsplice")
  data <- read.table(file_counts, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  counts <- data[, 2:7]
  tbl_exons <- table(sapply(strsplit(data$exon, ":"), function(s) s[[1]]))
  gene_IDs <- names(tbl_exons)
  n_exons <- unname(tbl_exons)
  condition <- rep(c("untreated", "treated"), each = 3)
  
  # run wrapper function
  # note: single core required for Travis CI
  # note: suppress warnings for grouped = FALSE in cv.glmnet due to small number of observations
  suppressWarnings(
    res <- regsplice(counts, gene_IDs, n_exons, condition, n_cores_reg = 1)
  )
  
  n_genes <- length(res@gene_IDs)
  
  n <- 81  # expect 81 genes remaining (out of 100 initially)
  
  
  expect_equal(n_genes, n)
  
  expect_length(res@p_vals, n)
  expect_length(res@p_adj, n)
  expect_length(res@LR_stats, n)
  expect_length(res@df_tests, n)
  
  expect_true(all(!is.na(res@p_vals)))
  expect_true(all(!is.na(res@p_adj)))
  
  expect_true(all(res@p_vals >= 0))
  expect_true(all(res@p_vals <= 1))
  expect_true(all(res@p_adj >= 0))
  expect_true(all(res@p_adj <= 1))
  
  expect_true(all(res@LR_stats >= 0, na.rm = TRUE))
  expect_true(all(res@df_tests >= 1, na.rm = TRUE))
})



