#' @include class_RegspliceData.R class_RegspliceResults.R
NULL




#' Fit models.
#' 
#' Model fitting functions for \code{regsplice} package.
#' 
#' There are three model fitting functions:
#' 
#' \code{fitRegMultiple} fits regularized (lasso) models containing an optimal subset of 
#' exon:condition interaction terms for each gene. The model fitting procedure penalizes 
#' the interaction terms only, so that the main effect terms for exons and samples are 
#' always included. This ensures that the null model is nested, allowing likelihood ratio
#' tests to be calculated.
#' 
#' \code{fitNullMultiple} fits the null models, which do not contain any interaction 
#' terms.
#' 
#' \code{fitFullMultiple} fits full models, which contain all exon:condition interaction
#' terms for each gene.
#' 
#' See \code{\link{createDesignMatrix}} for more details about the terms in each model.
#' 
#' The fitting functions fit models for all genes in the data set. The functions are 
#' parallelized using \code{BiocParallel::bplapply} for faster runtime. For 
#' \code{fitRegMultiple}, the default number of processor cores is 8, or the maximum 
#' available if less than 8. For \code{fitNullMultiple} and \code{fitFullMultiple}, the
#' default is one core, since these functions are already extremely fast for most data
#' sets.
#' 
#' A random seed can be provided with the \code{seed} argument, to generate reproducible 
#' results.
#' 
#' If the \code{rs_data} object does not contain a weights matrix, all exon bins are 
#' weighted equally.
#' 
#' Previous step: Initialize \code{\linkS4class{RegspliceResults}} object with 
#' \code{\link{initializeResults}}.
#' Next step: Calculate likelihood ratio tests with \code{\link{LRTests}}.
#' 
#' 
#' @param rs_results \code{\linkS4class{RegspliceResults}} object, which will be used to 
#'   store results. Initialized using the constructor function \code{RegspliceResults()}.
#'   See \code{\linkS4class{RegspliceResults}} for details.
#' @param rs_data \code{\linkS4class{RegspliceData}} object. In the case of RNA-seq read 
#'   count data, this has been pre-transformed with \code{\link{runVoom}}. Contains 
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
#' @param n_cores Number of cores for parallel evaluation. For \code{fitRegMultiple}, the
#'   default is 8, or the maximum available if less than 8. For \code{fitFullMultiple} 
#'   and \code{fitNullMultiple}, the default is 1, since these functions are already very
#'   fast.
#' @param seed Random seed (integer). Default is NULL. Provide an integer value to set 
#'   the random seed for reproducible results.
#' @param progress_bar Whether to display progress bar. Default is TRUE.
#' @param ... Other arguments to pass to \code{cv.glmnet}, \code{glmnet}, or \code{glm}.
#' 
#'   
#' @return Returns a \code{\linkS4class{RegspliceResults}} object containing the fitted 
#'   model objects, deviance of fitted models, and degrees of freedom of fitted models. 
#'   See \code{\linkS4class{RegspliceResults}} for details.
#' 
#' 
#' @seealso \code{\link{createDesignMatrix}} \code{\link{RegspliceResults}} 
#'   \code{\link{initializeResults}} \code{\link{LRTests}}
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
#' rs_data <- RegspliceData(counts, gene_IDs, n_exons, condition)
#' 
#' rs_data <- filterZeros(rs_data)
#' rs_data <- filterLowCounts(rs_data)
#' rs_data <- runNormalization(rs_data)
#' rs_data <- runVoom(rs_data)
#' 
#' rs_results <- initializeResults(rs_data)
#' 
#' rs_results <- fitRegMultiple(rs_results, rs_data, n_cores = 1)
#' rs_results <- fitNullMultiple(rs_results, rs_data)
#' rs_results <- fitFullMultiple(rs_results, rs_data)
#' 
fitRegMultiple <- function(rs_results, rs_data, 
                           alpha = 1, lambda_choice = c("lambda.min", "lambda.1se"), 
                           n_cores = NULL, seed = NULL, progress_bar = TRUE, ...) {
  
  lambda_choice <- match.arg(lambda_choice)
  
  gene_IDs <- names(table(rowData(rs_data)$gene_IDs))
  n_genes <- length(gene_IDs)
  
  FUN <- function(i) {
    gene_ID_i <- gene_IDs[i]
    if (gene_ID_i != rs_results@gene_IDs[i]) stop("gene IDs do not match")
    
    rs_data_i <- suppressMessages(rs_data[gene_ID_i, ])
    
    .fitRegSingle(rs_data = rs_data_i, alpha = alpha, lambda_choice = lambda_choice, ...)
  }
  
  message("Fitting regularized (lasso) models...")
  if (is.null(n_cores)) n_cores <- min(BiocParallel::multicoreWorkers(), 8)
  BPPARAM <- BiocParallel::MulticoreParam(workers = n_cores, 
                                          RNGseed = seed, 
                                          progressbar = progress_bar)
  
  # setting seed with BiocParallel when using glmnet doesn't work if using only one core;
  # use set.seed() instead
  if (n_cores == 1 & !is.null(seed)) set.seed(seed)
  
  out <- BiocParallel::bplapply(seq_len(n_genes), FUN = FUN, BPPARAM = BPPARAM)
  
  # collapse lists
  fit_collapse <- lapply(out, "[[", "fit")
  dev_collapse <- sapply(out, "[[", "dev")
  df_collapse <- sapply(out, "[[", "df")
  
  rs_results@fit_reg_models <- fit_collapse
  rs_results@fit_reg_dev <- dev_collapse
  rs_results@fit_reg_df <- df_collapse
  
  rs_results
}



#' @rdname fitRegMultiple
#' @export
#' 
fitNullMultiple <- function(rs_results, rs_data, n_cores = 1, 
                            seed = NULL, progress_bar = TRUE, ...) {
  
  gene_IDs <- names(table(rowData(rs_data)$gene_IDs))
  n_genes <- length(gene_IDs)
  
  FUN <- function(i) {
    gene_ID_i <- gene_IDs[i]
    if (gene_ID_i != rs_results@gene_IDs[i]) stop("gene IDs do not match")
    
    rs_data_i <- suppressMessages(rs_data[gene_ID_i, ])
    
    .fitNullSingle(rs_data_i, ...)
  }
  
  message("Fitting null models...")
  BPPARAM <- BiocParallel::MulticoreParam(workers = n_cores, 
                                          RNGseed = seed, 
                                          progressbar = progress_bar)
  
  out <- BiocParallel::bplapply(seq_len(n_genes), FUN = FUN, BPPARAM = BPPARAM)
  
  # collapse lists
  fit_collapse <- lapply(out, "[[", "fit")
  dev_collapse <- sapply(out, "[[", "dev")
  df_collapse <- sapply(out, "[[", "df")
  
  rs_results@fit_null_models <- fit_collapse
  rs_results@fit_null_dev <- dev_collapse
  rs_results@fit_null_df <- df_collapse
  
  rs_results
}



#' @rdname fitRegMultiple
#' @export
#' 
fitFullMultiple <- function(rs_results, rs_data, n_cores = 1, 
                            seed = NULL, progress_bar = TRUE, ...) {
  
  gene_IDs <- names(table(rowData(rs_data)$gene_IDs))
  n_genes <- length(gene_IDs)
  
  FUN <- function(i) {
    gene_ID_i <- gene_IDs[i]
    if (gene_ID_i != rs_results@gene_IDs[i]) stop("gene IDs do not match")
    
    rs_data_i <- suppressMessages(rs_data[gene_ID_i, ])
    
    .fitFullSingle(rs_data_i, ...)
  }
  
  message("Fitting full models...")
  BPPARAM <- BiocParallel::MulticoreParam(workers = n_cores, 
                                          RNGseed = seed, 
                                          progressbar = progress_bar)
  
  out <- BiocParallel::bplapply(seq_len(n_genes), FUN = FUN, BPPARAM = BPPARAM)
  
  # collapse lists
  fit_collapse <- lapply(out, "[[", "fit")
  dev_collapse <- sapply(out, "[[", "dev")
  df_collapse <- sapply(out, "[[", "df")
  
  rs_results@fit_full_models <- fit_collapse
  rs_results@fit_full_dev <- dev_collapse
  rs_results@fit_full_df <- df_collapse
  
  rs_results
}



