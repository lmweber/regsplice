library(regsplice)
context("Likelihood ratio tests function checks inputs")

test_that('fitted_models_GLM object is required if when_null_selected = "GLM"', {
  condition <- rep(c(0, 1), each = 3)
  n_exons <- 10
  Y <- list(as.data.frame(matrix(sample(100:200, 60, replace = TRUE), nrow = 10)))
  fitted_models_reg <- fit_reg(Y, condition)
  fitted_models_GLM <- fit_GLM(Y, condition)
  fitted_models_null <- fit_null(Y, condition)
  
  expect_error(LR_tests(fitted_models_reg = fitted_models_reg, 
                        fitted_models_null = fitted_models_null, 
                        when_null_selected = "GLM"))
})

