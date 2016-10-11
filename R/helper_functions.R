########################
### Helper functions ###
########################


# Remove genes containing only a single exon bin (since differential splicing requires
# multiple exon bins).
# 
.remove_single_exon_genes <- function(data) {
  
  if (!("RegspliceData" %in% is(data))) stop("'data' must be a 'RegspliceData' object")
  
  gene_IDs_rep <- rowData(data)$gene_IDs
  gene_IDs <- names(table(gene_IDs_rep))
  
  ix_singles <- rep(FALSE, length(gene_IDs))
  
  for (i in seq_along(gene_IDs)) {
    data_gene <- suppressMessages(data[gene_IDs[i], ])
    if (nrow(data_gene) == 1) {
      ix_singles[i] <- TRUE
    }
  }
  
  message(paste("removed", sum(ix_singles), "remaining single-exon gene(s)"))
  
  suppressMessages(data[gene_IDs[!ix_singles], ])
}


