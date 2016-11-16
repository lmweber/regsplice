library(regsplice)
context("voom weights and transformation/normalization")

test_that("voom weights and transformation/normalization are calculated correctly", {
  # prepared and filtered data for one gene
  Y <- list(data.frame(sample1 = c(1, 4, 0, 5), 
                       sample2 = c(0, 1, 0, 5), 
                       sample3 = c(0, 9, 6, 30), 
                       sample4 = c(3, 12, 4, 23), 
                       sample5 = c(1, 4, 0, 10), 
                       sample6 = c(1, 3, 1, 22)))
  names(Y) <- "gene12"
  condition <- rep(c(0, 1), each = 3)
  
  out_voom <- voom_weights(Y, condition)
  
  
  expect_length(out_voom, 2)
  
  expect_is(out_voom$Y[[1]], "data.frame")
  expect_is(out_voom$Y[[1]][2, 3], "numeric")
  
  expect_is(out_voom$weights[[1]], "data.frame")
  expect_is(out_voom$weights[[1]][2, 3], "numeric")
})

