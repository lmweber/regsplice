library(regsplice)
context("Vignette example")

test_that("results from vignette example are as expected", {
  
  # load data
  file_counts <- system.file("extdata/vignette_counts.txt", package = "regsplice")
  data <- read.table(file_counts, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  counts <- data[, 2:7]
  gene <- sapply(strsplit(data$exon, ":"), function(s) s[[1]])
  
  # create meta-data for biological samples
  condition <- rep(c("untreated", "treated"), each = 3)
  
  # prepare data
  Y <- prepare_data(counts = counts, gene = gene)
  
  # filter low-count exons
  Y <- filter_exons(Y = Y, filter_min_per_exon = 6, filter_min_per_sample = 3)
  
  # normalization factors
  norm_factors <- run_normalization(Y = Y)
  
  # voom transformation and weights
  out_voom <- run_voom(Y = Y, condition = condition, norm_factors = norm_factors)
  Y <- out_voom$Y
  weights <- out_voom$weights
  
  # fit models
  # note: single core required for Travis CI
  # note: suppress warnings for grouped = FALSE in cv.glmnet due to small number of observations
  suppressWarnings(
    fit_reg <- fit_models_reg(Y = Y, condition = condition, weights = weights, n_cores = 1)
  )
  fit_null <- fit_models_null(Y = Y, condition = condition, weights = weights)
  fit_GLM <- fit_models_GLM(Y = Y, condition = condition, weights = weights)
  
  # calculate likelihood ratio tests
  res <- LR_tests(fit_reg = fit_reg, 
                  fit_null = fit_null, 
                  when_null_selected = "ones")
  
  
  n_genes <- 81  # 81 genes after data preparation and filtering (out of 100)
  
  expect_length(Y, n_genes)
  expect_length(fit_reg$dev, n_genes)
  expect_length(fit_null$dev, n_genes)
  expect_length(fit_GLM$dev, n_genes)
  expect_length(res$p_val, n_genes)
  expect_length(res$p_adj, n_genes)
  expect_length(res$LR_stats, n_genes)
  expect_length(res$df_tests, n_genes)

  expect_true(all(!is.na(res$p_val)))
  expect_true(all(!is.na(res$p_adj)))

  expect_true(all(res$p_val >= 0))
  expect_true(all(res$p_val <= 1))
  expect_true(all(res$p_adj >= 0))
  expect_true(all(res$p_adj <= 1))

  expect_true(all(res$LR_stats >= 0, na.rm = TRUE))
  expect_true(all(res$df_tests >= 1, na.rm = TRUE))
})



test_that("results from vignette example are as expected (using wrapper function)", {
  
  # load data
  file_counts <- system.file("extdata/vignette_counts.txt", package = "regsplice")
  data <- read.table(file_counts, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  counts <- data[, 2:7]
  gene <- sapply(strsplit(data$exon, ":"), function(s) s[[1]])
  
  # create meta-data for biological samples
  condition <- rep(c("untreated", "treated"), each = 3)
  
  # run wrapper function
  # note: single core required for Travis CI
  # note: suppress warnings for grouped = FALSE in cv.glmnet due to small number of observations
  suppressWarnings(
    res <- regsplice(counts = counts, gene = gene, condition = condition, n_cores_reg = 1)
  )
  
  
  n_genes <- 81  # 81 genes after data preparation and filtering (out of 100)
  
  expect_length(res$p_val, n_genes)
  expect_length(res$p_adj, n_genes)
  expect_length(res$LR_stats, n_genes)
  expect_length(res$df_tests, n_genes)
  
  expect_true(all(!is.na(res$p_val)))
  expect_true(all(!is.na(res$p_adj)))
  
  expect_true(all(res$p_val >= 0))
  expect_true(all(res$p_val <= 1))
  expect_true(all(res$p_adj >= 0))
  expect_true(all(res$p_adj <= 1))
  
  expect_true(all(res$LR_stats >= 0, na.rm = TRUE))
  expect_true(all(res$df_tests >= 1, na.rm = TRUE))
})


