library(regsplice)
context("Vignette example")

test_that("results from vignette example are as expected", {
  
  # load data
  file_counts <- system.file("extdata/counts.txt", package = "regsplice")
  data <- read.table(file_counts, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  counts <- data[, 2:7]
  gene <- sapply(strsplit(data$exon, ":"), function(s) s[[1]])
  
  # create meta-data for biological samples
  condition <- rep(c("untreated", "treated"), each = 3)
  
  # prepare data
  Y <- prepare_data(counts = counts, gene = gene)
  
  # fit models
  fitted_models_reg <- fit_reg(Y = Y, condition = condition, n_cores = 1)
  fitted_models_GLM <- fit_GLM(Y = Y, condition = condition)
  fitted_models_null <- fit_null(Y = Y, condition = condition)
  
  # calculate likelihood ratio tests
  res <- LR_tests(fitted_models_reg = fitted_models_reg, 
                  fitted_models_GLM = NULL, 
                  fitted_models_null = fitted_models_null, 
                  when_null_selected = "ones")
  
  
  n_genes <- 87  # 87 genes after data preparation and filtering (out of 100)
  
  expect_length(Y, n_genes)
  expect_length(fitted_models_reg$dev, n_genes)
  expect_length(fitted_models_GLM$dev, n_genes)
  expect_length(fitted_models_null$dev, n_genes)
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



test_that("results from vignette example are as expected (using wrapper function)", {
  
  # load data
  file_counts <- system.file("extdata/counts.txt", package = "regsplice")
  data <- read.table(file_counts, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  counts <- data[, 2:7]
  gene <- sapply(strsplit(data$exon, ":"), function(s) s[[1]])
  
  # create meta-data for biological samples
  condition <- rep(c("untreated", "treated"), each = 3)
  
  # run wrapper function
  res <- regsplice(counts = counts, gene = gene, condition = condition, n_cores_reg = 1)
  
  
  n_genes <- 87  # 87 genes after data preparation and filtering (out of 100)
  
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


