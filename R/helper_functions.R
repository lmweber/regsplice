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
  gene_IDs_singles <- gene_IDs[table(gene_IDs_rep) == 1]
  
  ix_singles <- gene_IDs_rep %in% gene_IDs_singles
  
  message(paste("removed", sum(ix_singles), "remaining single-exon gene(s)"))
  
  suppressMessages(rs_data[!ix_singles, ])
}


