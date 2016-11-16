#' @include class_RegspliceData.R class_RegspliceResults.R
NULL




#' Fit models.
#' 
#' Model fitting functions for \code{regsplice} package.
#' 
#' There are three model fitting functions:
#' 
#' \code{fit_reg_multiple} fits regularized (lasso) models containing an optimal subset
#' of exon:condition interaction terms for each gene. The model fitting procedure
#' penalizes the interaction terms only, so that the main effect terms for exons and
#' samples are always included. This ensures that the null model is nested, allowing
#' likelihood ratio tests to be calculated.
#' 
#' \code{fit_null_multiple} fits the null models, which do not contain any interaction 
#' terms.
#' 
#' \code{fit_full_multiple} fits full models, which contain all exon:condition
#' interaction terms for each gene.
#' 
#' See \code{\link{create_design_matrix}} for more details about the terms in each model.
#' 
#' The fitting functions fit models for all genes in the data set. The functions are 
#' parallelized using \code{BiocParallel::bplapply} for faster runtime. For 
#' \code{fit_reg_multiple}, the default number of processor cores is 8, or the maximum 
#' available if less than 8. For \code{fit_null_multiple} and \code{fit_full_multiple},
#' the default is one core, since these functions are already extremely fast for most
#' data sets.
#' 
#' A random seed can be provided with the \code{seed} argument, to generate reproducible
#' results.
#' 
#' If the \code{data} object does not contain a weights matrix, all exon bins are
#' weighted equally.
#' 
#' Previous step: Initialize \code{\linkS4class{RegspliceResults}} object with
#' \code{\link{initialize_results}}.
#' Next step: Calculate likelihood ratio tests with \code{\link{LR_tests}}.
#' 
#' 
#' @param results \code{\linkS4class{RegspliceResults}} object, which will be used to
#'   store results. Initialized using the constructor function \code{RegspliceResults()}.
#'   See \code{\linkS4class{RegspliceResults}} for details.
#' @param data \code{\linkS4class{RegspliceData}} object. In the case of RNA-seq read 
#'   count data, this has been pre-transformed with \code{\link{run_voom}}. Contains 
#'   \code{counts} and \code{weights} data matrices, and vector of experimental 
#'   conditions for each biological sample stored in \code{colData}. See 
#'   \code{\linkS4class{RegspliceData}} for details.
#' @param alpha Elastic net parameter \code{alpha} for \code{glmnet} model fitting 
#'   functions. Must be between 0 (ridge regression) and 1 (lasso). Default is 1 (lasso).
#'   See \code{glmnet} documentation for more details.
#' @param lambda_choice Parameter to select which optimal \code{lambda} value to choose 
#'   from the \code{cv.glmnet} cross validation fit. Choices are "lambda.min" (model with
#'   minimum cross-validated error) and "lambda.1se" (most regularized model with 
#'   cross-validated error within one standard error of minimum). Default is 
#'   "lambda.min". See \code{glmnet} documentation for more details.
#' @param n_cores Number of cores for parallel evaluation. For \code{fit_reg_multiple},
#'   the default is 8, or the maximum available if less than 8. For
#'   \code{fit_full_multiple} and \code{fit_null_multiple}, the default is 1, since these
#'   functions are already very fast.
#' @param seed Random seed (integer). Default is NULL. Provide an integer value to set 
#'   the random seed for reproducible results.
#' @param progress_bar Whether to display progress bar (\code{fit_reg_multiple} only). 
#'   Default is TRUE.
#' @param ... Other arguments to pass to \code{cv.glmnet}, \code{glmnet}, or \code{glm}.
#' 
#'   
#' @return Returns a \code{\linkS4class{RegspliceResults}} object containing the fitted 
#'   model objects, deviance of fitted models, and degrees of freedom of fitted models.
#'   See \code{\linkS4class{RegspliceResults}} for details.
#' 
#' 
#' @seealso \code{\link{create_design_matrix}} \code{\link{RegspliceResults}}
#'   \code{\link{initialize_results}} \code{\link{LR_tests}}
#' 
#' @seealso \code{\link[glmnet]{glmnet}} \code{\link[glmnet]{cv.glmnet}} 
#'   \code{\link[stats]{glm}}
#'   
#' @importFrom BiocParallel multicoreWorkers MulticoreParam bplapply
#' 
#' @export
#' 
#' @examples
#' file_counts <- system.file("extdata/vignette_counts.txt", package = "regsplice")
#' data <- read.table(file_counts, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#' head(data)
#' 
#' counts <- data[, 2:7]
#' tbl_exons <- table(sapply(strsplit(data$exon, ":"), function(s) s[[1]]))
#' gene_IDs <- names(tbl_exons)
#' n_exons <- unname(tbl_exons)
#' condition <- rep(c("untreated", "treated"), each = 3)
#' 
#' Y <- RegspliceData(counts, gene_IDs, n_exons, condition)
#' 
#' Y <- filter_zeros(Y)
#' Y <- filter_low_counts(Y)
#' Y <- run_normalization(Y)
#' Y <- run_voom(Y)
#' 
#' res <- initialize_results(Y)
#' 
#' res <- fit_reg_multiple(res, Y, n_cores = 1)
#' res <- fit_null_multiple(res, Y)
#' res <- fit_full_multiple(res, Y)
#' 
fit_reg_multiple <- function(results, data, 
                           alpha = 1, lambda_choice = c("lambda.min", "lambda.1se"), 
                           n_cores = NULL, seed = NULL, progress_bar = TRUE, ...) {
  
  lambda_choice <- match.arg(lambda_choice)
  
  gene_IDs <- names(table(rowData(data)$gene_IDs))
  n_genes <- length(gene_IDs)
  
  FUN <- function(i) {
    gene_ID_i <- gene_IDs[i]
    if (gene_ID_i != results@gene_IDs[i]) stop("gene IDs do not match")
    
    data_i <- suppressMessages(data[gene_ID_i, ])
    
    fit_reg_single(data = data_i, alpha = alpha, lambda_choice = lambda_choice, ...)
  }
  
  message("Fitting regularized (lasso) models...")
  if (is.null(n_cores)) n_cores <- min(BiocParallel::multicoreWorkers(), 8)
  BPPARAM <- BiocParallel::MulticoreParam(workers = n_cores, 
                                          RNGseed = seed, 
                                          progressbar = progress_bar)
  
  # setting seed with BiocParallel when using glmnet doesn't work if using only one core;
  # use set.seed() instead
  if (n_cores == 1 & !is.null(seed)) set.seed(seed)
  
  res <- BiocParallel::bplapply(seq_len(n_genes), FUN = FUN, BPPARAM = BPPARAM)
  
  # collapse lists
  fit_collapse <- lapply(res, "[[", "fit")
  dev_collapse <- sapply(res, "[[", "dev")
  df_collapse <- sapply(res, "[[", "df")
  
  results@fit_reg_models <- fit_collapse
  results@fit_reg_dev <- dev_collapse
  results@fit_reg_df <- df_collapse
  
  results
}



#' @rdname fit_reg_multiple
#' @export
#' 
fit_null_multiple <- function(results, data, n_cores = 1, seed = NULL, ...) {
  
  gene_IDs <- names(table(rowData(data)$gene_IDs))
  n_genes <- length(gene_IDs)
  
  FUN <- function(i) {
    gene_ID_i <- gene_IDs[i]
    if (gene_ID_i != results@gene_IDs[i]) stop("gene IDs do not match")
    
    data_i <- suppressMessages(data[gene_ID_i, ])
    
    fit_null_single(data_i, ...)
  }
  
  message("Fitting null models...")
  BPPARAM <- BiocParallel::MulticoreParam(workers = n_cores, RNGseed = seed)
  
  res <- BiocParallel::bplapply(seq_len(n_genes), FUN = FUN, BPPARAM = BPPARAM)
  
  # collapse lists
  fit_collapse <- lapply(res, "[[", "fit")
  dev_collapse <- sapply(res, "[[", "dev")
  df_collapse <- sapply(res, "[[", "df")
  
  results@fit_null_models <- fit_collapse
  results@fit_null_dev <- dev_collapse
  results@fit_null_df <- df_collapse
  
  results
}



#' @rdname fit_reg_multiple
#' @export
#' 
fit_full_multiple <- function(results, data, n_cores = 1, seed = NULL, ...) {
  
  gene_IDs <- names(table(rowData(data)$gene_IDs))
  n_genes <- length(gene_IDs)
  
  FUN <- function(i) {
    gene_ID_i <- gene_IDs[i]
    if (gene_ID_i != results@gene_IDs[i]) stop("gene IDs do not match")
    
    data_i <- suppressMessages(data[gene_ID_i, ])
    
    fit_full_single(data_i, ...)
  }
  
  message("Fitting full models...")
  BPPARAM <- BiocParallel::MulticoreParam(workers = n_cores, RNGseed = seed)
  
  res <- BiocParallel::bplapply(seq_len(n_genes), FUN = FUN, BPPARAM = BPPARAM)
  
  # collapse lists
  fit_collapse <- lapply(res, "[[", "fit")
  dev_collapse <- sapply(res, "[[", "dev")
  df_collapse <- sapply(res, "[[", "df")
  
  results@fit_full_models <- fit_collapse
  results@fit_full_dev <- dev_collapse
  results@fit_full_df <- df_collapse
  
  results
}



