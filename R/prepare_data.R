# Function to prepare data

# split up data into list of lists, where each list item contains data for one gene


prepare_data <- function(counts, gene) {
  # use data frame so each sub-matrix keeps its shape
  if (is.matrix(counts)) counts <- as.data.frame(counts)
  if (is.null(colnames(counts))) colnames(counts) <- paste0("sample", 1:ncol(counts))
  
  # keep genes (levels) in same order as provided in input
  if (!is.character(gene)) gene <- as.character(gene)
  gene <- factor(gene, levels = unique(gene))
  
  # split data by unique gene identifiers
  Y <- split(counts, gene)
  
  # remove any genes with zero expression
  zeros <- sapply(Y, function(d) all(d == 0))
  Y <- Y[!zeros]
}

