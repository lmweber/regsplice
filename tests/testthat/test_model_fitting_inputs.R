library(regsplice)
context("Model fitting functions check inputs")

test_that("providing data for a single gene as a vector returns an error", {
  condition <- rep(c(0, 1), each = 3)
  Y <- sample(100:200, 60, replace = TRUE)
  
  expect_error(fit_reg_single(Y = Y, condition = condition))
  expect_error(fit_null_single(Y = Y, condition = condition))
  expect_error(fit_GLM_single(Y = Y, condition = condition))
})


test_that("providing data for multiple genes as a single data frame returns an error", {
  condition <- rep(c(0, 1), each = 3)
  Y <- as.data.frame(matrix(sample(100:200, 60, replace = TRUE), nrow = 10))
  
  expect_error(fit_models_reg(Y = Y, condition = condition))
  expect_error(fit_models_null(Y = Y, condition = condition))
  expect_error(fit_models_GLM(Y = Y, condition = condition))
})

