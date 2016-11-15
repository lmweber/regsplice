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
#' The fitting functions fit models for all genes in the data set.
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
#' @param seed Random seed (integer). Default is NULL. Provide an integer value to set 
#'   the random seed for reproducible results.
#' @param ... Other arguments to pass to \code{cv.glmnet}, \code{glmnet}, or \code{glm}.
#' 
#'   
#' @return Returns a \code{\linkS4class{RegspliceResults}} object containing deviance and
#'   degrees of freedom of the fitted models. See \code{\linkS4class{RegspliceResults}}
#'   for details.
#' 
#' 
#' @seealso \code{\link{createDesignMatrix}} \code{\link{RegspliceResults}} 
#'   \code{\link{initializeResults}} \code{\link{LRTests}}
#' 
#' @seealso \code{\link[glmnet]{glmnet}} \code{\link[glmnet]{cv.glmnet}} 
#'   \code{\link[stats]{glm}}
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
#' rs_results <- fitRegMultiple(rs_results, rs_data)
#' rs_results <- fitNullMultiple(rs_results, rs_data)
#' rs_results <- fitFullMultiple(rs_results, rs_data)
#' 
fitRegMultiple <- function(rs_results, rs_data, 
                           alpha = 1, lambda_choice = c("lambda.min", "lambda.1se"), 
                           seed = NULL, ...) {
  
  lambda_choice <- match.arg(lambda_choice)
  
  gene_IDs <- names(table(rowData(rs_data)$gene_IDs))
  n_genes <- length(gene_IDs)
  
  FUN <- function(i) {
    gene_ID_i <- gene_IDs[i]
    if (gene_ID_i != rs_results@gene_IDs[i]) stop("gene IDs do not match")
    
    rs_data_i <- rs_data[gene_ID_i, ]
    
    .fitRegSingle(rs_data = rs_data_i, alpha = alpha, lambda_choice = lambda_choice, ...)
  }
  
  message("Fitting regularized (lasso) models...")
  
  if (is.null(seed)) seed <- as.numeric(Sys.time())
  set.seed(seed)
  
  out <- lapply(seq_len(n_genes), FUN = FUN)
  
  # collapse lists
  dev_collapse <- sapply(out, "[[", "dev")
  df_collapse <- sapply(out, "[[", "df")
  
  rs_results@fit_reg_dev <- dev_collapse
  rs_results@fit_reg_df <- df_collapse
  
  rs_results
}



#' @rdname fitRegMultiple
#' @export
#' 
fitNullMultiple <- function(rs_results, rs_data, seed = NULL, ...) {
  
  gene_IDs <- names(table(rowData(rs_data)$gene_IDs))
  n_genes <- length(gene_IDs)
  
  FUN <- function(i) {
    gene_ID_i <- gene_IDs[i]
    if (gene_ID_i != rs_results@gene_IDs[i]) stop("gene IDs do not match")
    
    rs_data_i <- rs_data[gene_ID_i, ]
    
    .fitNullSingle(rs_data_i, ...)
  }
  
  message("Fitting null models...")
  
  if (is.null(seed)) seed <- as.numeric(Sys.time())
  set.seed(seed)
  
  out <- lapply(seq_len(n_genes), FUN = FUN)
  
  # collapse lists
  dev_collapse <- sapply(out, "[[", "dev")
  df_collapse <- sapply(out, "[[", "df")
  
  rs_results@fit_null_dev <- dev_collapse
  rs_results@fit_null_df <- df_collapse
  
  rs_results
}



#' @rdname fitRegMultiple
#' @export
#' 
fitFullMultiple <- function(rs_results, rs_data, seed = NULL, ...) {
  
  gene_IDs <- names(table(rowData(rs_data)$gene_IDs))
  n_genes <- length(gene_IDs)
  
  FUN <- function(i) {
    gene_ID_i <- gene_IDs[i]
    if (gene_ID_i != rs_results@gene_IDs[i]) stop("gene IDs do not match")
    
    rs_data_i <- rs_data[gene_ID_i, ]
    
    .fitFullSingle(rs_data_i, ...)
  }
  
  message("Fitting full models...")
  
  if (is.null(seed)) seed <- as.numeric(Sys.time())
  set.seed(seed)
  
  out <- lapply(seq_len(n_genes), FUN = FUN)
  
  # collapse lists
  dev_collapse <- sapply(out, "[[", "dev")
  df_collapse <- sapply(out, "[[", "df")
  
  rs_results@fit_full_dev <- dev_collapse
  rs_results@fit_full_df <- df_collapse
  
  rs_results
}



