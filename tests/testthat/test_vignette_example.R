library(regsplice)
context("Vignette example")

test_that("example from vignette gives correct p-values", {
  
  # load data
  file_counts <- system.file("extdata/counts.txt", package = "regsplice")
  data <- read.table(file_counts, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  counts <- data[, 2:7]
  gene <- sapply(strsplit(data$exon, ":"), function(s) s[[1]])
  
  # create meta-data for biological samples
  condition <- rep(c("untreated", "treated"), each = 3)
  
  # prepare data
  Y <- split_genes(counts = counts, gene = gene)
  Y <- filter_genes(Y)
  
  # fit models (note random seed required)
  fitted_models_reg <- fit_reg(Y = Y, condition = condition, seed = 123)
  fitted_models_GLM <- fit_GLM(Y = Y, condition = condition)
  fitted_models_null <- fit_null(Y = Y, condition = condition)
  
  # calculate likelihood ratio tests
  res <- LR_tests(fitted_models_reg = fitted_models_reg, 
                  fitted_models_GLM = NULL, 
                  fitted_models_null = fitted_models_null, 
                  when_null_selected = "ones")
  
  # saved p-values
  file_saved <- system.file("tests/testthat/p_vals_vignette.txt", package = "regsplice")
  p_vals_saved <- read.table(file_saved, header = TRUE)
  p_vals_saved <- unname(unlist(p_vals_saved))
  
  expect_equal(res$p_vals, p_vals_saved)
})

