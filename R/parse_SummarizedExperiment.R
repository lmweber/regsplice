# Internal function to process SummarizedExperiment object passed as input to
# regsplice() wrapper function.
# 
parse_SummarizedExperiment <- function(se) {
  
  if (!is(se, "SummarizedExperiment")) {
    stop("input is not a SummarizedExperiment")
  }
  
  if (length(assays(se)) > 1) {
    warning("SummarizedExperiment input object contains more than one 'assays' object. ", 
            "Only the first 'assays' object will be used; this is assumed to contain the read counts.")
  }
  
  counts <- assays(se)$counts
  gene_IDs <- rowData(se)$gene_IDs
  n_exons <- unname(table(rowData(se)$gene_IDs))
  condition <- colData(se)$condition
  
  if (!all(gene_IDs == names(table(rowData(se)$gene_IDs)))) {
    stop("gene IDs do not match (this may be due to alphabetical re-ordering)")
  }
  
  if (is.null(counts)) {
    stop("'counts' matrix could not be identified in SummarizedExperiment input object")
  }
  
  if (is.null(gene_IDs)) {
    stop("'gene_IDs' vector could not be found in rowData of SummarizedExperiment input object")
  }
  
  if (is.null(n_exons)) {
    stop("number of exons per gene could not be calculated from repeated entries in 'gene_IDs' vector")
  }
  
  if (is.null(condition)) {
    stop("'condition' vector could not be found in colData of SummarizedExperiment input object")
  }
  
  RegspliceData(counts, gene_IDs, n_exons, condition)
}


