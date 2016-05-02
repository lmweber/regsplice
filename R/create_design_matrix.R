#####################################
## Functions for regsplice package ##
## author: Lukas Weber             ##
#####################################

# Create design matrix for a single gene.



#' Create design matrix.
#' 
#' Creates a gene-level design matrix in the format required by the \code{regsplice}
#' package.
#' 
#' Creates a design matrix for a single gene, given the number of exons and group 
#' membership for each sample. This function will typically be called once for every gene
#' in a data set.
#' 
#' Notes:
#' 
#' \itemize{
#'   \item No main effect is included for the group factor. The group main effect is not 
#'   required since it is absorbed within the sample main effects.
#'   \item No intercept is included, since it is simpler to let the \code{glmnet}
#'   functions add it back later.
#' }
#' 
#' @param group Vector containing group identifiers for each sample.
#' @param nexons Number of exons in the gene.
#'   
#' @return Returns a design matrix in the format required by the functions in the
#'   \code{regsplice} package.
#'   
#' @examples
#' group <- rep(c(0, 1), each = 3)
#' nexons <- 10
#' createDesignMatrix(group, nexons)
create_design_matrix <- function(condition, n_exons) {
  n_samples <- length(condition)
  
  Exon <- factor(rep(1:n_exons, times = n_samples))
  Samp <- factor(rep(1:n_samples, each = n_exons))
  Cond <- factor(rep(condition, each = n_exons))
  
  # Build design matrix manually, since model.matrix(~ Exon + Samp + Exon:Cond)[, -1] 
  # includes an extra (linearly dependent) interaction column for the first exon.
  # Also don't include an intercept column, since easier if model fitting functions add 
  # it back later.
  main_effects <- model.matrix(~ Exon + Samp)[, -1, drop = FALSE]
  int_temp     <- model.matrix(~ Exon + Cond + Exon:Cond)[, -1, drop = FALSE]
  interactions <- int_temp[, grep(":", colnames(int_temp)), drop = FALSE]
  
  X <- cbind(main_effects, interactions)
  
  X
}


