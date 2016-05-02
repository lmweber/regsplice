library(regsplice)
context("Vignette example")

test_that("vignette example gives correct results", {
  # load data
  file_counts <- system.file("extdata/counts.txt", package = "regsplice")
  data <- read.table(file_counts, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  counts <- data[, 2:7]
  gene <- sapply(strsplit(data$exon, ":"), function(s) s[[1]])
  
  # prepare data
  Y <- prepare_data(counts, gene)
  
  # fit models
  fit_reg  <- fit_models_reg(Y = Y, condition = condition)  ### can use n_cores = 4 on Mac
  fit_GLM  <- fit_models_GLM(Y = Y, condition = condition)
  fit_null <- fit_models_null(Y = Y, condition = condition)
  
  # calculate likelihood ratio tests
  LR_tests(fit_reg = fit_reg, fit_GLM = fit_GLM, fit_null = fit_null, when_null_selected = "ones")
  
  
})