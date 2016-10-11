library(regsplice)
context("Filter low-count exons")

test_that("low-count exon bins are filtered correctly", {
  
  counts <- rbind(c(1, 0, 0, 0, 0, 0), 
                  c(2, 3, 5, 1, 2, 1), 
                  c(0, 2, 1, 1, 2, 2), 
                  c(0, 4, 0, 0, 0, 0), 
                  c(0, 3, 4, 1, 1, 2), 
                  c(10, 5, 6, 8, 2, 8), 
                  c(13, 8, 9, 2, 4, 2))
  gene_IDs <- paste0("gene", 1:3)
  n_exons <- c(3, 2, 2)
  condition <- rep(c(0, 1), each = 3)
  
  rs_data <- RegspliceData(counts, gene_IDs, n_exons, condition)
  
  rs_data <- filterZeros(rs_data)
  rs_data <- filterLowCounts(rs_data)
  
  n_genes <- length(names(table(rowData(rs_data)$gene_IDs)))
  n_exons <- nrow(rowData(rs_data))
  
  expect_equivalent(n_genes, 1)
  expect_equivalent(n_exons, 2)
})



