#####################################
## Functions for regsplice package ##
## author: Lukas Weber             ##
#####################################


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
createDesignMatrix <- function(group, nexons) {
  nsamples <- length(group)
  
  Exon <- factor( rep(1:nexons, nsamples) )
  Grp  <- factor( rep(group, each=nexons) )
  Samp <- factor( rep(1:nsamples, each=nexons) )
  
  # build matrix by first including Grp main effect and then removing it - otherwise 
  # model.matrix keeps an extra (linearly dependent) interaction column
  X <- model.matrix(~ Exon + Samp + Grp + Exon:Grp)[, -1]
  X <- X[, !(colnames(X) %in% paste0("Grp", levels(Grp)))]
  
  return(X)
}

