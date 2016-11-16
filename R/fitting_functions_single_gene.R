# Model fitting functions for a single gene. These are internal functions, which should 
# not be called directly by users.


#' @importFrom glmnet cv.glmnet deviance.glmnet
#' @importFrom methods is
#' 
fit_reg_single <- function(data, alpha = 1, 
                           lambda_choice = c("lambda.min", "lambda.1se"), ...) {
  
  if (!("RegspliceData" %in% is(data))) stop("'data' must be a 'RegspliceData' object")
  
  Y <- countsData(data)
  weights <- weightsData(data)
  condition <- colData(data)$condition
  
  n_exons <- nrow(Y)
  X <- create_design_matrix(condition = condition, n_exons = n_exons)
  if (is.null(weights)) weights <- matrix(1, nrow = n_exons, ncol = length(condition))
  
  lambda_choice <- match.arg(lambda_choice)
  
  # convert Y and weights to vectors (note matrix is read by column)
  Y <- as.vector(Y)
  weights <- as.vector(weights)
  
  # identify interaction columns by ":" in column names
  int_cols <- grepl(":", colnames(X))
  
  # which columns to penalize during lasso fitting
  pen <- as.numeric(int_cols)
  
  # Note: 'while' loop and 'tryCatch' statement are included to repeat the model fitting
  # procedure if 'cv.glmnet' returns an error, due to the bug introduced in 'glmnet'
  # version 2.0-2. See unit test 'test_glmnet_error_example.R' or GitHub repository at 
  # https://github.com/lmweber/glmnet-error-example for details. Since the error depends 
  # on the random seed, we simply re-run the model fitting procedure (with changing 
  # random seed) until it passes (maximum of 20 tries).
  fit <- NULL
  i <- 0
  while (is.null(fit) & i < 20) {
    tryCatch(fit <- glmnet::cv.glmnet(x = X, y = Y, weights = weights, 
                                      alpha = alpha, penalty.factor = pen, 
                                      intercept = TRUE, standardize = FALSE, ...), 
             error = function(e) {
               # print any unrelated error messages, which do not contain the term 'predmat'
               if (!grepl("predmat", e)) print(as.character(e))
             })
    i <- i + 1
  }
  
  if (is.null(fit)) {
    dev <- NA
    df <- NA
  } else {
    ix_opt <- which(fit$lambda == fit[[lambda_choice]])
    dev <- glmnet::deviance.glmnet(fit$glmnet.fit)[ix_opt]
    df <- fit$glmnet.fit$df[ix_opt]
  }
  
  list(dev = dev, df = df, fit = fit)
}



#' @importFrom stats glm
#' @importFrom methods is
#' 
fit_null_single <- function(data, ...) {
  
  if (!("RegspliceData" %in% is(data))) stop("'data' must be a 'RegspliceData' object")
  
  Y <- countsData(data)
  weights <- weightsData(data)
  condition <- colData(data)$condition
  
  n_exons <- nrow(Y)
  X <- create_design_matrix(condition = condition, n_exons = n_exons)
  if (is.null(weights)) weights <- matrix(1, nrow = n_exons, ncol = length(condition))
  
  # identify interaction columns by ":" in column names
  int_cols <- grepl(":", colnames(X))
  
  # convert Y and weights to vectors (note matrix is read by column)
  Y <- as.vector(Y)
  weights <- as.vector(weights)
  
  fit <- stats::glm(Y ~ X[, !int_cols], weights = weights, ...)
  
  dev <- fit$deviance
  df <- fit$df.null - fit$df.residual
  
  list(dev = dev, df = df, fit = fit)
}



#' @importFrom stats glm
#' @importFrom methods is
#' 
fit_full_single <- function(data, ...) {
  
  if (!("RegspliceData" %in% is(data))) stop("'data' must be a 'RegspliceData' object")
  
  Y <- countsData(data)
  weights <- weightsData(data)
  condition <- colData(data)$condition
  
  n_exons <- nrow(Y)
  X <- create_design_matrix(condition = condition, n_exons = n_exons)
  if (is.null(weights)) weights <- matrix(1, nrow = n_exons, ncol = length(condition))
  
  # convert Y and weights to vectors (note matrix is read by column)
  Y <- as.vector(Y)
  weights <- as.vector(weights)
  
  fit <- stats::glm(Y ~ X, weights = weights, ...)
  
  dev <- fit$deviance
  df <- fit$df.null - fit$df.residual
  
  list(dev = dev, df = df, fit = fit)
}


