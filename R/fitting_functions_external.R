#####################################
## Functions for regsplice package ##
## author: Lukas Weber             ##
#####################################


#' Fit models for multiple genes.
#' 
#' Parallelized functions to fit models for multiple genes.
#' 
#' \code{fitRegModel} fits regularized models, \code{fitGLM} fits full GLMs containing 
#' all interaction terms, and \code{fitNull} fits null models with no interaction terms. 
#' These functions are parallelized to take advantage of multiple processors or cores.
#' 
#' @section Regularized models:
#' \code{fitRegModel} fits l1-regularized (lasso) models for each gene using the 
#' \code{glmnet} package. Together with the null models from \code{fitNullModel}, these
#' can then be used to calculate likelihood ratio p-values using the function 
#' \code{\link{lrTest}}.
#' 
#' The l1-regularization (lasso) model fitting procedure in this function is designed to 
#' penalize interaction terms only. The fitted models therefore consist of all main
#' effect terms for exons and samples, together with some or all of the interaction
#' terms. The nested null models consist of the main effect terms only. Likelihood ratio
#' tests can therefore be calculated to compare the evidence for the two models.
#' 
#' For some genes, the regularized model will typically select zero interaction terms. In
#' these cases the fitted and null models are equivalent, so a likelihood ratio test
#' cannot be calculated. The \code{\link{lrTest}} function provides two possible
#' strategies for these cases — setting p-values for these genes to 1, or re-fitting a
#' full GLM containing all interaction terms. See documentation for \code{\link{lrTest}}
#' for more details.
#' 
#' Notes:
#' 
#' \itemize{
#'   Note that interaction terms are identified by the column names of the design matrices 
#' in \code{X_genes}. Column names containing a colon \code{:} are assumed to represent 
#' interaction columns (for example \code{Exon2:Grp1}). Design matrices in this format
#' can be generated with the function \code{\link{createDesignMatrix}}. Alternatively,
#' the arguments \code{group} and \code{nexons_genes} can be provided instead of
#' \code{X_genes} (see below).
#' 
#'   The expression data matrices in \code{Y_genes} are assumed to be on a continuous
#' scale. RNA-seq count data can be transformed to a continuous scale using the
#' \code{voom} function in the \code{limma} package.
#' 
#'   The call to \code{glmnet} sets the argument \code{standardize = FALSE}. 
#' Standardization is not required here since the design matrices in \code{X_genes} 
#' contain only indicator variables.
#' }
#' 
#' 
#' @section Full models containing all interaction terms:
#' \code{fitGLM} fits full GLMs containing all interaction terms for each gene. The model
#' fitting is done with \code{glm}. These models can also be used together with the null
#' models from \code{\link{fitNullModels}} to calculate likelihood ratio tests with
#' \code{\link{lrTest}}.
#' 
#' Alternatively, the full models from \code{fitGLM} can be used to replace fitted
#' regularized models containing zero interaction terms during the calculation of the
#' likelihood ratio tests wiht \code{\link{lrTests}}. This option can be selected with
#' the argument \code{when_null_selected = "GLM"} in the call to \code{\link{lrTest}} 
#' (see documentation for \code{\link{lrTest}} for more details).
#' 
#' The full models contain exon main effects, sample main effects, and all exon:group 
#' interaction terms. This is similar to the approach used in \code{voom-diffSlice} 
#' (function \code{diffSplice} in the \code{limma} package).
#' 
#' 
#' @section Null models:
#' \code{fitNullModel} fits null models with no interaction terms for each gene. The 
#' model fitting is done with \code{glm}. The null models can be used together with the
#' fitted regularized or full models to calculate likelihood ratio tests with the
#' function \code{\link{lrTest}}.
#' 
#' The null models contain exon main effects and sample main effects, but do not contain
#' any interaction terms. Hence they are nested within the regularized models from
#' \code{\link{fitRegModel}} or the full models from \code{\link{GLM}}, allowing
#' likelihood ratio tests to be calculated.
#' 
#' 
#' @section Parallelization:
#' To speed up runtime by taking advantage of multiple processors or cores, these
#' functions are are parallelized using the \code{MulticoreParam} function from the
#' \code{BiocParallel} package. This works with Mac OSX and Linux systems — on Windows
#' systems, the code will simply run on a single core. The parallelization is implemented
#' in the internal function \code{\link{fitParallel}}.
#' 
#' 
#' @param X_genes List of design matrices. Individual design matrices in the required 
#'   format can be created with the function \code{\link{createDesignMatrix}}. Note that 
#'   column names are used to identify interaction terms (see above). Can be set to
#'   \code{NULL}, in which case \code{group} and \code{nexons_genes} must be provided.
#' @param Y_genes List of expression data matrices. Assumed to be on a continuous scale, 
#'   for example transformed using the \code{voom} function from the \code{limma}
#'   package.
#' @param weights_genes Optional list of weights matrices, for example generated by
#'   \code{voom}.
#' @param group Optional vector of group membership identifiers for each sample. Must be
#'   provided if \code{X_genes = NULL}, so that \code{\link{createDesignMatrix}} can
#'   be called to generate design matrices.
#' @param nexons_genes Optional vector containing number of exons for each gene. Must be
#'   provided if \code{X_genes = NULL}, so that \code{\link{createDesignMatrix}} can
#'   be called to generate design matrices.
#' @param alpha Elastic net parameter \code{alpha} for \code{glmnet} fitting functions.
#'   Required for \code{fitRegModel}. Must be between 0 (ridge regression) and 1 (lasso).
#'   See \code{glmnet} documentation for details. Default is 1 (lasso).
#' @param lambda_choice Parameter to select which optimal lambda value to choose from the
#'   \code{cv.glmnet} cross validation fit. Required for \code{fitRegModel}. Choices are
#'   "lambda.min" (model with minimum cross-validated error) and "lambda.1se" (most
#'   regularized model with cross-validated error within one standard error of minimum).
#'   See \code{glmnet} documentation for more details. Default is "lambda.min".
#' @param return_fitted Logical; whether to return fitted model objects from
#'   \code{cv.glmnet}.
#' @param ncores Number of cores for parallel evaluation. Default is 1. Note that on
#'   Windows systems, the code will always run on a single core.
#' @param ... Other arguments to pass to \code{cv.glmnet}, \code{glmnet}, or \code{glm}.
#' 
#' @return Returns a list containing:
#' \itemize{
#'   fit_genes: (list) (optional) fitted model objects (\code{fitGLM} or
#'   \code{fitNullModel}) or optimal fitted model objects (\code{fitRegModel}) for each
#'   gene
#'   dev_genes: (vector) deviance of fitted models (\code{fitGLM} or \code{fitNullModel})
#'   or optimal fitted models (\code{fitRegModel}) for each gene
#'   df_genes: (vector) degrees of freedom of fitted models (\code{fitGLM} or
#'   \code{fitNullModel}) or optimal fitted models (\code{fitRegModel}) for each gene
#' }
#' 
#' @seealso \code{\link{lrTest}} \code{\link[glmnet]{glmnet}} 
#'   \code{\link[glmnet]{cv.glmnet}} \code{\link[stats]{glm}}
#' 
#' @examples
#' set.seed(1)
#' group <- rep(c(0, 1), each=3)
#' nexons <- 8
#' X <- createDesignMatrix(group, nexons)
#' Y <- rnorm(nrow(X), mean=2, sd=1)
#' ix <- c(7, 8) + ( 8 * rep(0:5, each=2) )
#' Y[ix] <- Y[ix] + 1
#' fitRegModel(X, Y)
#' fitGLM(X, Y)
#' fitNullModel(X, Y)
#' 
fitRegModel <- function(X_genes = NULL, Y_genes, weights_genes = NULL, 
                        group = NULL, nexons_genes = NULL, 
                        alpha = 1, lambda_choice = c("lambda.min", "lambda.1se"), 
                        return_fitted = FALSE, ncores = 1, ...) {
  
  fitParallel(fitting_function = "regularized", 
              X_genes = X_genes, Y_genes = Y_genes, weights_genes = weights_genes, 
              group = group, nexons_genes = nexons_genes, 
              alpha = alpha, lambda_choice = lambda_choice, 
              return_fitted = return_fitted, ncores = ncores, ...)
}



#' @rdname fitRegModel
#' 
fitGLM <- function(X_genes = NULL, Y_genes, weights_genes = NULL, 
                   group = NULL, nexons_genes = NULL, 
                   return_fitted = FALSE, ncores = 1, ...) {
  
  fitParallel(fitting_function = "GLM", 
              X_genes = X_genes, Y_genes = Y_genes, weights_genes = weights_genes, 
              group = group, nexons_genes = nexons_genes, 
              return_fitted = return_fitted, ncores = ncores, ...)
}



#' @rdname fitRegModel
#' 
fitNullModel <- function(X_genes = NULL, Y_genes, weights_genes = NULL, 
                         group = NULL, nexons_genes = NULL, 
                         return_fitted = FALSE, ncores = 1, ...) {
  
  fitParallel(fitting_function = "null", 
              X_genes = X_genes, Y_genes = Y_genes, weights_genes = weights_genes, 
              group = group, nexons_genes = nexons_genes, 
              return_fitted = return_fitted, ncores = ncores, ...)
}


