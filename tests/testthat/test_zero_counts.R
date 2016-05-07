library(regsplice)
context("Genes with zero counts are removed")

test_that("genes with zero counts are correctly removed", {
  counts <- rbind(matrix(1, nrow = 5, ncol = 4), 
                  matrix(0, nrow = 3, ncol = 4), 
                  matrix(1, nrow = 2, ncol = 4))
  gene <- paste0("gene", rep(1:4, times = c(3, 2, 3, 2)))
  
  Y <- split_genes(counts = counts, gene = gene)
  Y_filtered <- filter_genes(Y = Y)
  
  expect_length(Y, 4)
  expect_length(Y_filtered, 3)
})

