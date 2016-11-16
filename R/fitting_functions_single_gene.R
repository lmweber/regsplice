# Model fitting functions for a single gene. These are internal functions, which should 
# not be called directly by users.


#' @importFrom glmnet cv.glmnet deviance.glmnet
#' 
fit_reg_single <- function(Y, condition, weights = NULL, alpha = 1, 
                           lambda_choice = c("lambda.min", "lambda.1se"), ...) {
  
  if (!(is.matrix(Y) | is.data.frame(Y))) {
    stop("data Y for a single gene must be a matrix or data frame")
  }
  if (is.data.frame(Y)) {
    Y <- as.matrix(Y)
    row.names(Y) <- 1:nrow(Y)
  }
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
  # random seed) until it passes.
  fit <- NULL
  while (is.null(fit)) {
    tryCatch(fit <- glmnet::cv.glmnet(x = X, y = Y, weights = weights, 
                                      alpha = alpha, penalty.factor = pen, 
                                      intercept = TRUE, standardize = FALSE, ...), 
             error = function(e) {
               # print any unrelated error messages (which do not contain the term 'predmat')
               if (!grepl("predmat", e)) print(as.character(e))
             })
  }
  
  ix_opt <- which(fit$lambda == fit[[lambda_choice]])
  dev <- glmnet::deviance.glmnet(fit$glmnet.fit)[ix_opt]
  df <- fit$glmnet.fit$df[ix_opt]
  
  list(dev = dev, df = df, fit = fit)
}



#' @importFrom stats glm
#' 
fit_GLM_single <- function(Y, condition, weights = NULL, ...) {
  
  if (!(is.matrix(Y) | is.data.frame(Y))) {
    stop("data Y for a single gene must be a matrix or data frame")
  }
  if (is.data.frame(Y)) {
    Y <- as.matrix(Y)
    row.names(Y) <- 1:nrow(Y)
  }
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



#' @importFrom stats glm
#' 
fit_null_single <- function(Y, condition, weights = NULL, ...) {
  
  if (!(is.matrix(Y) | is.data.frame(Y))) {
    stop("data Y for a single gene must be a matrix or data frame")
  }
  if (is.data.frame(Y)) {
    Y <- as.matrix(Y)
    row.names(Y) <- 1:nrow(Y)
  }
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


