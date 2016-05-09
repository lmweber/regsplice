library(regsplice)
context("Design matrices")

test_that("design matrix is generated correctly for saved example", {
  condition <- rep(c(0, 1), each = 3)
  n_exons <- 4
  design <- create_design_matrix(condition = condition, n_exons = n_exons)
  
  file_saved <- system.file("extdata/design_matrix_example.txt", package = "regsplice")
  design_saved <- as.matrix(read.table(file_saved, header = TRUE, check.names = FALSE))
  
  expect_equivalent(design, design_saved)
})

