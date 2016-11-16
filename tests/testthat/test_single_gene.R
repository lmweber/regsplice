library(regsplice)
context("Test fitting functions for a single gene work")

test_that("testing", {
  condition <- rep(c(0, 1), each = 3)
  Y <- matrix(sample(20:49, 60, replace = TRUE), nrow = 10)
  
  fit_reg_model_single(Y = Y, condition = condition)
  fit_GLM_single(Y = Y, condition = condition)
  fit_null_model_single(Y = Y, condition = condition)
})