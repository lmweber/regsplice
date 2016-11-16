# Function to prepare data

# create SummarizedExperiment object


prepare_data <- function(counts, gene, exon = NULL, condition) {
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  if (!is.character(gene)) gene <- as.character(gene)
  if (!is.null(exon) & !is.character(exon)) exon <- as.character(exon)
  if (!is.character(condition)) condition <- as.character(condition)
  
  if (is.null(colnames(counts))) colnames(counts) <- paste0("sample", 1:ncol(counts))
  
  se <- SummarizedExperiment::SummarizedExperiment(assays = counts)
  
  # add meta-data
  
  if (is.null(exon)) {
    row_mdata <- S4Vectors::DataFrame(gene = gene)
  } else {
    row_mdata <- S4Vectors::DataFrame(gene = gene, exon = exon)
  }
  S4Vectors::mcols(se) <- row_mdata
  
  col_mdata <- S4Vectors::DataFrame(condition = condition)
  SummarizedExperiment::colData(se) <- col_mdata
  
  se
}

