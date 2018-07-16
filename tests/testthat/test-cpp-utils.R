context("util.cpp unit tests")

# Test c++ function add_constant -----------------------------------------------
test_that("add_constant", {
  v1 <- c(1, 2, 3)
  v2 <- c(1.2, 4, 5.6)
  expect_equal(hesim:::C_test_add_constant_int(v1, 5.5),
               floor(v1 + 5.5))
  expect_equal(hesim:::C_test_add_constant_double(v2, 5.5),
               v2 + 5.5)
})
