library(regsplice)
context("glmnet error example")

test_that("glmnet error example passes", {
  
  # 'glmnet' version 2.0-2 introduced a bug in the cross validation step, which gives 
  # errors in the 'regsplice' lasso model fitting function '.fitRegSingle()' for some 
  # random seeds. Earlier versions of 'glmnet' (<= v. 1.9-8) did not have this problem. 
  # The 'glmnet' package authors have advised that the bug will be fixed in a future 
  # update. Until then, we have included checks in '.fitRegSingle()' to simply re-run 
  # the model fitting procedure if the error occurs. This unit test checks that the 
  # 'regsplice' functions pass for a specific example of the bug. For more details, see 
  # the GitHub repository at: https://github.com/lmweber/glmnet-error-example
  
  # LOAD DATA
  # Y: vector of response values
  # X: design matrix of indicator variables and interaction terms
  file <- system.file("extdata/data_glmnet_error_example.txt", package = "regsplice")
  data <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE)
  Y <- data[, 1]
  X <- as.matrix(data[, -1])
  
  # PENALTY FACTOR
  # penalize interaction terms only
  penalty <- rep(c(0, 1), times = c(13, 8))
  
  # SHOW ERROR BY RUNNING CV.GLMNET DIRECTLY (not run)
  # cv.glmnet: cross validation for lambda
  # The following code demonstrates the error.
  # error message: 'Error in predmat[which, seq(nlami)] = preds : replacement has length zero'
  #set.seed(1)
  #glmnet::cv.glmnet(x = X, y = Y, family = "gaussian", penalty.factor = penalty)
  
  # RUN REGSPLICE MODEL FITTING FUNCTION
  # regsplice model fitting function checks for the error and should pass
  counts <- matrix(Y, ncol = 6)
  gene_IDs <- "gene1"
  n_exons <- 9
  condition <- rep(c(0, 1), each = 3)
  
  rs_data <- RegspliceData(counts, gene_IDs, n_exons, condition)
  
  # filtering low-count exons can also solve the problem
  #rs_data <- filterZeros(rs_data)
  #rs_data <- filterLowCounts(rs_data)
  
  rs_results <- initializeResults(rs_data)
  
  # error occurs for set.seed(1) if checks are not included in '.fitRegSingle()'
  set.seed(1)
  rs_results <- fitRegMultiple(rs_results, rs_data)
  
  
  # no output object is returned if the error occurs
  expect_is(rs_results, "RegspliceResults")
  expect_is(rs_results@fit_reg_models, "list")
  expect_length(rs_results@fit_reg_models, 1)
  
})



