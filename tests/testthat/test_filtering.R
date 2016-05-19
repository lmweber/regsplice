library(regsplice)
context("Data filtering")

test_that("low-count exons are filtered correctly", {
  counts <- rbind(c(1, 0, 0, 0, 0, 0), 
                  c(2, 3, 5, 1, 2, 1), 
                  c(0, 2, 1, 1, 2, 2), 
                  c(0, 4, 0, 0, 0, 0), 
                  c(0, 3, 4, 1, 1, 2), 
                  c(10, 5, 6, 8, 2, 8))
  gene <- paste0("gene", rep(1:2, times = c(3, 3)))
  
  Y <- prepare_data(counts = counts, gene = gene)
  Y <- filter_exons(Y = Y)
  
  n_exons <- sapply(Y, nrow)
  
  expect_length(Y, 1)
  expect_equivalent(n_exons, 2)
})

