#####################################
## Functions for regsplice package ##
## author: Lukas Weber             ##
#####################################


#' Fit model for a single gene.
#' 
#' Internal functions to fit models for a single gene.
#' 
#' \code{fitRegModelSingle} fits a regularized model.
#' 
#' \code{fitGLMSingle} fits a full GLM containing all interaction terms.
#' 
#' \code{fitNullModelSingle} fits a null model with no interaction terms.
#' 
#' Each function fits a model for a single gene only. These are intended as internal
#' functions. In most cases, the user will only need to access the external functions
#' \code{\link{fitRegModel}}, \code{\link{fitGLM}}, and \code{\link{fitNullModel}}, which
#' can fit models for multiple genes using parallelization.
#' 
#' The external functions also call the internal function \code{\link{fitParallel}}.
#' 
#' See documentation for \code{\link{fitRegModel}}, \code{\link{fitGLM}}, and
#' \code{\link{fitNullModel}} for more details, including parameters and return values.
#' 
#' @seealso \code{\link{fitRegModel}} \code{\link{fitGLM}} \code{\link{fitNullModel}} 
#'   \code{\link{fitParallel}} \code{\link{lrTest}} \code{\link[glmnet]{glmnet}} 
#'   \code{\link[glmnet]{cv.glmnet}} \code{\link[stats]{glm}}
#'   
fit_reg_model_single <- function(Y, condition, weights = NULL, alpha = 1, 
                                 lambda_choice = c("lambda.min", "lambda.1se"), ...) {
  
  if (!is.matrix(Y)) stop("data Y for a single gene must be provided as a matrix")
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



#' @rdname fitRegModelSingle
#' 
fit_GLM_single <- function(Y, condition, weights = NULL, ...) {
  
  if (!is.matrix(Y)) stop("data Y for a single gene must be provided as a matrix")
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



#' @rdname fitRegModelSingle
#' 
fit_null_model_single <- function(Y, condition, weights = NULL, ...) {
  
  if (!is.matrix(Y)) stop("data Y for a single gene must be provided as a matrix")
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


