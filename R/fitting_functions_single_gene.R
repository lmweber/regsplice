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
fitRegModelSingle <- function(X = NULL, Y, weights = NULL, 
                              group = NULL, nexons = NULL, 
                              alpha = 1, lambda_choice = c("lambda.min", "lambda.1se"), 
                              ...) {
  
  if (is.null(X)) X <- createDesignMatrix(group = group, nexons = nexons)
  lambda_choice <- match.arg(lambda_choice)
  if (is.null(weights)) weights <- rep(1, length(Y))
  
  # identify interaction columns
  int_cols <- grepl(":", colnames(X))
  
  fit <- glmnet::cv.glmnet(x = X, y = Y, weights = weights, 
                           alpha = alpha, penalty.factor = as.numeric(int_cols), 
                           standardize = FALSE, ...)
  
  ix_opt <- which(fit$lambda == fit[[lambda_choice]])
  dev <- deviance(fit$glmnet.fit)[ix_opt]
  df <- fit$glmnet.fit$df[ix_opt]
  
  return(list(fit = fit, dev = dev, df = df))
}



#' @rdname fitRegModelSingle
#' 
fitGLMSingle <- function(X = NULL, Y, weights = NULL, 
                         group = NULL, nexons = NULL, ...) {
  
  if (is.null(X)) X <- createDesignMatrix(group = group, nexons = nexons)
  if (is.null(weights)) weights <- rep(1, length(Y))
  
  # identify interaction columns
  int_cols <- grepl(":", colnames(X))
  
  fit <- glm(Y ~ X, weights = weights, ...)
  dev <- fit$deviance
  df <- fit$df.null - fit$df.residual
  
  return(list(fit = fit, dev = dev, df = df))
}



#' @rdname fitRegModelSingle
#' 
fitNullModelSingle <- function(X = NULL, Y, weights = NULL, 
                               group = NULL, nexons = NULL, ...) {
  
  if (is.null(X)) X <- createDesignMatrix(group = group, nexons = nexons)
  if (is.null(weights)) weights <- rep(1, length(Y))
  
  # identify interaction columns
  int_cols <- grepl(":", colnames(X))
  
  fit <- glm(Y ~ X[, !int_cols], weights = weights, ...)
  dev <- fit$deviance
  df <- fit$df.null - fit$df.residual
  
  return(list(fit = fit, dev = dev, df = df))
}


