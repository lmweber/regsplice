#' Create design matrix.
#' 
#' Create a model design matrix for a single gene.
#' 
#' Creates a model design matrix for a single gene in the format required by the 
#' \emph{regsplice} model fitting functions. Required inputs are the biological
#' conditions for each sample, and the number of exons in the gene.
#' 
#' The design matrix includes main effect terms for each exon and each sample, and
#' interaction terms between the exons and the biological conditions.
#' 
#' Note that the design matrix does not include main effect terms for the biological
#' conditions, since these are absorbed into the main effect terms for the samples. In
#' addition, the design matrix does not include an intercept column, since it is simpler
#' to let the model fitting functions add an intercept term later.
#' 
#' The model fitting functions in subsequent steps will call this function once for each
#' gene.
#' 
#' @param condition Vector containing the biological conditions for each sample 
#'   (character or numeric vector, or factor).
#' @param n_exons Number of exons in the gene (integer).
#' 
#' @return Returns a model design matrix for the gene, in the format required by the 
#'   \emph{regsplice} model fitting functions.
#' 
#' @family create_design_matrix fit_reg fit_GLM fit_null LR_tests
#' 
#' @export
#'   
#' @examples
#' condition <- rep(c(0, 1), each = 3)
#' n_exons <- 10
#' create_design_matrix(condition, n_exons)
#' 
create_design_matrix <- function(condition, n_exons) {
  
  n_samples <- length(condition)
  
  Exon <- factor(rep(1:n_exons, times = n_samples))
  Samp <- factor(rep(1:n_samples, each = n_exons))
  Cond <- factor(rep(condition, each = n_exons))
  
  # Build design matrix manually, since creating it automatically with 
  # model.matrix(~ Exon + Samp + Exon:Cond)[, -1] will include an extra (linearly
  # dependent) interaction column for the first exon. Also don't include an intercept
  # column, since it is simpler to let the model fitting functions add it back later.
  main_effects <- stats::model.matrix(~ Exon + Samp)[, -1, drop = FALSE]
  int_temp     <- stats::model.matrix(~ Exon + Cond + Exon:Cond)[, -1, drop = FALSE]
  interactions <- int_temp[, grep(":", colnames(int_temp)), drop = FALSE]
  
  X <- cbind(main_effects, interactions)
  
  X
}

