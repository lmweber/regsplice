########################
### Helper functions ###
########################


# Remove genes containing only a single exon bin (since differential splicing requires
# multiple exon bins).
# 
.removeSingleExonGenes <- function(rs_data) {
  
  if (!("RegspliceData" %in% is(rs_data))) stop("'rs_data' must be a 'RegspliceData' object")
  
  gene_IDs_rep <- rowData(rs_data)$gene_IDs
  gene_IDs <- names(table(gene_IDs_rep))
  
  ix_singles <- rep(FALSE, length(gene_IDs))
  
  for (i in seq_along(gene_IDs)) {
    rs_data_gene <- suppressMessages(rs_data[gene_IDs[i], ])
    if (nrow(rs_data_gene) == 1) {
      ix_singles[i] <- TRUE
    }
  }
  
  message(paste("removed", sum(ix_singles), "remaining single-exon gene(s)"))
  
  suppressMessages(rs_data[gene_IDs[!ix_singles], ])
}


