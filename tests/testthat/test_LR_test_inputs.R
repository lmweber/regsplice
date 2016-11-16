library(regsplice)
context("Likelihood ratio tests function checks inputs")

test_that('fit_GLM object is required if when_null_selected = "GLM"', {
  condition <- rep(c(0, 1), each = 3)
  n_exons <- 10
  Y <- list(as.data.frame(matrix(sample(100:200, 60, replace = TRUE), nrow = 10)))
  
  fit_reg  <- fit_models_reg(Y, condition)
  fit_null <- fit_models_null(Y, condition)
  fit_GLM  <- fit_models_GLM(Y, condition)
  
  expect_error(LR_tests(fit_reg = fit_reg, 
                        fit_null = fit_null, 
                        when_null_selected = "GLM"))
})

