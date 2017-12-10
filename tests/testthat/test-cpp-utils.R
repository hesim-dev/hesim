context("util.cpp unit tests")

# Test c++ function matrix_byrow ----------------------------------------------
test_that("matrix_byrow", {
  vec <- c(1, 2, 3, 4, 5, 6)
  matR <- matrix(vec, nrow = 2, ncol = 3, byrow = TRUE)
  matC <- matrix_byrow(vec, 2, 3)
  expect_equal(matR, matC)
})

# Test c++ function matrix_bycol ----------------------------------------------
test_that("matrix_bycol", {
  vec <- c(1, 2, 3, 4, 5, 6)
  matR <- matrix(vec, nrow = 2, ncol = 3, byrow = FALSE)
  matC <- matrix_bycol(vec, 2, 3)
  expect_equal(matR, matC)
})


