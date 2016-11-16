library(regsplice)
context("subsetting RegspliceData objects")

test_that("row-subsetting RegspliceData objects", {
  
  set.seed(123)
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  gene_IDs <- paste0("gene", 1:length(n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  rs_data <- RegspliceData(counts, gene_IDs, n_exons, condition)
  
  # row-subsetting using row index
  expect_equal(nrow(rs_data[1, ]), 1)
  expect_equal(nrow(rs_data[1:3, ]), 3)
  
  # row-subsetting using gene IDs
  expect_equal(nrow(rs_data["gene1", ]), 7)
  expect_equal(nrow(rs_data["gene2", ]), 18)
})




test_that("column-subsetting RegspliceData objects", {
  
  set.seed(123)
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  gene_IDs <- paste0("gene", 1:length(n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  rs_data <- RegspliceData(counts, gene_IDs, n_exons, condition)
  
  # column-subsetting
  expect_equal(ncol(rs_data[, 1:2]), 2)
  expect_equal(ncol(rs_data[, 4:6]), 3)
  expect_equal(ncol(rs_data[, 6]), 1)
})




test_that("row-and-column-subsetting RegspliceData objects", {
  
  set.seed(123)
  n_exons <- c(7, 18, 5, 5, 1, 3, 11)
  counts <- matrix(sample(100:200, sum(n_exons) * 6, replace = TRUE), ncol = 6)
  gene_IDs <- paste0("gene", 1:length(n_exons))
  condition <- rep(c(0, 1), each = 3)
  
  rs_data <- RegspliceData(counts, gene_IDs, n_exons, condition)
  
  # subsetting both rows and columns
  
  expect_equal(nrow(rs_data[1, 1]), 1)
  expect_equal(ncol(rs_data[1, 1]), 1)
  
  expect_equal(nrow(rs_data[2, 3]), 1)
  expect_equal(ncol(rs_data[2, 3]), 1)
  
  expect_equal(nrow(rs_data[5:10, 4:6]), 6)
  expect_equal(ncol(rs_data[5:10, 4:6]), 3)
  
  expect_equal(nrow(rs_data["gene1", 1:3]), 7)
  expect_equal(ncol(rs_data["gene1", 1:3]), 3)
  
  expect_equal(nrow(rs_data[c("gene3", "gene6"), 4]), 8)
  expect_equal(ncol(rs_data[c("gene3", "gene6"), 4]), 1)
  
  expect_equal(nrow(rs_data[c("gene1", "gene2"), 1:3]), 25)
  expect_equal(ncol(rs_data[c("gene1", "gene2"), 1:3]), 3)
})



