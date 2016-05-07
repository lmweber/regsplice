#' Create design matrix.
#' 
#' Creates a design matrix for a single gene.
#' 
#' Creates a design matrix for a single gene in the format required by the
#' \emph{regsplice} fitting functions. Required inputs are the conditions for each
#' biological sample, and the number of exons in the gene.
#' 
#' The fitting functions in subsequent steps will call this function once for each gene.
#' 
#' Notes:
#' 
#' \itemize{
#' \item The design matrix does not contain any main effect terms for the conditions, 
#' since these are absorbed within the main effect terms for the samples.
#' 
#' \item The design matrix does not include an intercept column, since it is simpler to 
#' let the model fitting functions add it back later.
#' }
#' 
#' @param condition Vector containing the biological conditions for each sample 
#'   (character or numeric vector, or factor).
#' @param n_exons Number of exons in the gene (integer).
#' 
#' @return Returns a design matrix for the gene, in the format required by the
#'   \emph{regsplice} model fitting functions.
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
  
  # Build design matrix manually, since model.matrix(~ Exon + Samp + Exon:Cond)[, -1] 
  # includes an extra (linearly dependent) interaction column for the first exon. Also 
  # don't include an intercept column, since it is simpler to let the model fitting 
  # functions add it back later.
  main_effects <- model.matrix(~ Exon + Samp)[, -1, drop = FALSE]
  int_temp     <- model.matrix(~ Exon + Cond + Exon:Cond)[, -1, drop = FALSE]
  interactions <- int_temp[, grep(":", colnames(int_temp)), drop = FALSE]
  
  X <- cbind(main_effects, interactions)
  
  X
}

