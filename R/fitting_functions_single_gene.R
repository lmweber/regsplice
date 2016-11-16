# Fitting functions for a single gene. These are internal functions, which are not
# exported.


fit_reg_model_single <- function(Y, condition, weights = NULL, alpha = 1, 
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
  pen <- as.numeric(int_cols)
  
  fit <- glmnet::cv.glmnet(x = X, y = Y, weights = weights, 
                           alpha = alpha, penalty.factor = pen, 
                           intercept = TRUE, standardize = FALSE, ...)
  
  ix_opt <- which(fit$lambda == fit[[lambda_choice]])
  dev <- deviance(fit$glmnet.fit)[ix_opt]
  df <- fit$glmnet.fit$df[ix_opt]
  
  list(fit = fit, dev = dev, df = df)
}



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
  
  fit <- glm(Y ~ X, weights = weights, ...)
  dev <- fit$deviance
  df <- fit$df.null - fit$df.residual
  
  list(fit = fit, dev = dev, df = df)
}



fit_null_model_single <- function(Y, condition, weights = NULL, ...) {
  
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
  
  fit <- glm(Y ~ X[, !int_cols], weights = weights, ...)
  dev <- fit$deviance
  df <- fit$df.null - fit$df.residual
  
  list(fit = fit, dev = dev, df = df)
}


