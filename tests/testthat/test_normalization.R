library(regsplice)
context("Normalization factors")

test_that("normalization factors work correctly", {
  
  counts <- rbind(c(101, 198, 53, 102, 100, 101), 
                  c(90, 220, 60, 89, 90, 103), 
                  c(95, 205, 66, 105, 85, 99), 
                  c(300, 622, 145, 311, 290, 320), 
                  c(290, 602, 140, 290, 303, 298), 
                  c(290, 410, 90, 190, 211, 209))
  gene <- paste0("gene", rep(1:2, times = c(3, 3)))
  
  Y <- prepare_data(counts = counts, gene = gene)
  Y <- filter_exons(Y = Y)
  
  norm_factors <- run_normalization(Y)
  
  expect_length(norm_factors, 6)
  expect_is(norm_factors, "numeric")
})


