###########################
### S4 class definition ###
###########################


#' @rdname RegspliceData
#' @export
#' 
setClass("RegspliceData", contains = "SummarizedExperiment")




############################
### Constructor function ###
############################


#' RegspliceData objects.
#' 
#' \code{RegspliceData} objects contain data in the format required by functions in the
#' \code{regsplice} package.
#' 
#' The \code{RegspliceData} format is based on the 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} container. Initially, objects
#' contain raw data along with meta-data for rows (genes and exons) and columns 
#' (biological samples). During subsequent steps in the \code{regsplice} analysis 
#' pipeline, the data values are modified, and additional data and meta-data are added to
#' the object. Final results are stored in a \code{\linkS4class{RegspliceResults}} 
#' object.
#' 
#' \code{RegspliceData} objects are created with the constructor function
#' \code{RegspliceData()}. Arguments for the constructor function are described above.
#' 
#' After creating the \code{RegspliceData} object, the function 
#' \code{\link{filter_zeros}} can be used to begin pre-processing the data.
#' 
#' 
#' @param counts RNA-seq read counts or exon microarray intensities (matrix or data 
#'   frame). Rows are exons, and columns are biological samples.
#' @param gene_IDs Vector of gene IDs (character vector). Length is equal to the number 
#'   of genes.
#' @param n_exons Vector of exon lengths (numeric vector of integers), i.e. the number of
#'   exon bins per gene. Length is equal to the number of genes.
#' @param condition Experimental condition for each biological sample (character or 
#'   numeric vector, or factor).
#' @param x \code{RegspliceData} object (for accessor or subsetting functions).
#' @param i Gene names (character vector) or row numbers (numeric vector) for subsetting 
#'   genes or exons. Note that when subsetting whole genes, gene names (character vector)
#'   should be provided instead of row numbers, to avoid possible errors due to selecting
#'   incorrect row numbers. Row numbers may be provided to subset individual exons.
#' @param j Column numbers (numeric vector) for subsetting biological samples.
#' @param ... Additional arguments for replacement with \code{`[<-`}.
#' @param value Value for replacement with \code{`[<-`}.
#' @param withDimnames See \code{SummarizedExperiment::assays()}.
#' 
#' 
#' @field counts Matrix of RNA-seq read counts or exon microarray intensities. Rows are 
#'   exons, and columns are biological samples.
#' @field weights (Optional) Matrix of observation-level weights. Rows are exons, and 
#'   columns are biological samples. Created by the \code{\link{run_voom}} function.
#' @field rowData \code{DataFrame} of row meta-data. This should contain two columns: 
#'   \code{gene_IDs} and \code{exon_IDs}, which are created by the \code{RegspliceData}
#'   constructor function.
#' @field colData \code{DataFrame} of column meta-data. This contains the experimental 
#'   condition and (optionally) normalization factors for each biological sample. 
#'   Normalization factors are created by the \code{\link{run_voom}} function.
#' 
#' 
#' @section Accessor functions:
#' 
#' \itemize{
#' \item \code{countsData()}: Accesses the \code{counts} data matrix.
#' \item \code{weightsData()}: Accesses the (optional) \code{weights} data matrix.
#' \item \code{rowData()}: Accesses the \code{DataFrame} of row meta-data. This should
#' contain two columns: \code{gene_IDs} and \code{exon_IDs}.
#' \item \code{colData()}: Accesses the \code{DataFrame} of column meta-data. This
#' contains the experimental condition and (optionally) normalization factors for each
#' biological sample.
#' }
#' 
#' 
#' @section Subsetting:
#' 
#' Subsetting of \code{RegspliceData} objects is performed with square brackets, 
#' \code{x[i, j]}, where \code{x} is the name of the object. The subsetting operations
#' are designed to keep data and meta-data in sync.
#' 
#' For subsetting by rows, there are two possibilities:
#' \itemize{
#' \item Subsetting genes: To subset whole genes, provide a character vector of gene 
#' names to the argument \code{i}. The returned object will contain all rows 
#' corresponding to these genes. Row numbers should not be used when subsetting whole
#' genes, since this risks potential errors due to selecting incorrect rows.
#' \item Subsetting exons: To subset individual exons, provide the corresponding row
#' numbers to the argument \code{i}.
#' }
#' 
#' For subsetting by columns (biological samples), provide the corresponding column
#' numbers to the argument \code{j}.
#' 
#' 
#' @return Returns a \code{RegspliceData} object.
#' 
#' @seealso \code{\link{filter_zeros}}
#' 
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom SummarizedExperiment SummarizedExperiment Assays
#' @importFrom methods new
#' 
#' @export
#' 
#' @examples
#' # ---------
#' # Example 1
#' # ---------
#' 
#' counts <- matrix(sample(100:200, 14 * 6, replace = TRUE), nrow = 14, ncol = 6)
#' gene_IDs <- paste0("gene", 1:5)
#' n_exons <- c(3, 2, 3, 1, 5)
#' condition <- rep(c(0, 1), each = 3)
#' 
#' Y <- RegspliceData(counts, gene_IDs, n_exons, condition)
#' 
#' Y
#' countsData(Y)
#' rowData(Y)
#' colData(Y)
#' 
#' Y[1, ]
#' Y[1, 1:3]
#' 
#' Y["gene1", ]
#' Y["gene1", 1:3]
#' 
#' 
#' # --------------------
#' # Example 2 (Vignette)
#' # --------------------
#' 
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
#' Y
#' head(countsData(Y))
#' rowData(Y)
#' colData(Y)
#' 
#' Y[1, ]
#' Y[1, 1:3]
#' 
#' Y["ENSG00000000003", ]
#' Y["ENSG00000000003", 1:3]
#' 
RegspliceData <- function(counts, gene_IDs, n_exons, condition) {
  
  if (!(is.matrix(counts) | is.data.frame(counts))) {
    stop("'counts' must be a matrix or data frame")
  }
  if (!is.character(gene_IDs)) {
    stop("'gene_IDs' must be a character vector")
  }
  if (!is.numeric(n_exons)) {
    stop("'n_exons' must be a numeric vector")
  }
  if (!(is.character(condition) | is.numeric(condition) | is.factor(condition))) {
    stop("'condition' must be a character vector, numeric vector, or factor")
  }
  
  if (sum(n_exons) != nrow(counts)) {
    stop("total number of exons 'sum(n_exons)' does not match number of rows in counts")
  }
  
  gene_IDs_rep <- unname(unlist(mapply(rep, gene_IDs, n_exons)))
  
  # generate exon IDs as 3-digit numbers with leading zeros
  exon_IDs <- sprintf("%03d", unlist(lapply(n_exons, seq)))
  
  rowData <- S4Vectors::DataFrame(gene_IDs = gene_IDs_rep, exon_IDs = exon_IDs)
  
  if (!is.null(colnames(counts))) {
    sample_names <- colnames(counts)
    colData <- S4Vectors::DataFrame(sample_names = sample_names, condition = condition)
  } else {
    colData <- S4Vectors::DataFrame(condition = condition)
  }
  
  counts <- SummarizedExperiment::Assays(S4Vectors::SimpleList(counts = as.matrix(counts)))
  
  new("RegspliceData", 
      assays = counts, 
      elementMetadata = rowData, 
      colData = colData)
}




##########################
### Accessor functions ###
##########################


#' @rdname RegspliceData
#' @importFrom SummarizedExperiment assays "assays<-"
#' @importFrom methods callNextMethod
#' @export
#' 
setMethod("assays", "RegspliceData", function(x, ..., value, withDimnames) {
  # note: "value" is required for setter operation "assays<-"
  callNextMethod(x, ..., value, withDimnames)
})




#' @rdname RegspliceData
#' @export
#' 
setGeneric("countsData", function(x) {
  standardGeneric("countsData")
})


#' @rdname RegspliceData
#' @importFrom SummarizedExperiment assays
#' @export
#' 
setMethod("countsData", "RegspliceData", function(x) {
  SummarizedExperiment::assays(x)$counts
})




#' @rdname RegspliceData
#' @export
#' 
setGeneric("weightsData", function(x) {
  standardGeneric("weightsData")
})


#' @rdname RegspliceData
#' @importFrom SummarizedExperiment assays
#' @export
#' 
setMethod("weightsData", "RegspliceData", function(x) {
  SummarizedExperiment::assays(x)$weights
})




#' @rdname RegspliceData
#' @importFrom SummarizedExperiment rowData
#' @importFrom methods callNextMethod
#' @export
#' 
setMethod("rowData", "RegspliceData", function(x) {
  callNextMethod(x)
})




#' @rdname RegspliceData
#' @importFrom SummarizedExperiment colData "colData<-"
#' @importFrom methods callNextMethod
#' @export
#' 
setMethod("colData", "RegspliceData", function(x, ..., value) {
  # note: "value" is required for setter operation "colData<-"
  callNextMethod(x, ..., value)
})




###########################
### Subsetting function ###
###########################


#' @rdname RegspliceData
#' @importFrom methods callNextMethod
#' @export
#' 
setMethod("[", "RegspliceData", function(x, i, j) {

  if (is.character(i)) {
    message("subsetting rows by gene names")
    
    keep_i <- rowData(x)$gene_IDs %in% i
    callNextMethod(x = x, i = keep_i, j = j)

  } else {
    message("subsetting rows by row numbers")
    
    callNextMethod()
  }
})



