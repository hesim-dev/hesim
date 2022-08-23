context("rtrunc.cpp unit tests")

test_that("rtrunc_repeat", {
  lower <- 5
  upper <- 8
  r <- replicate(10, hesim:::C_test_rtrunc_repeat(lower, upper))
  expect_true(all(r >= lower & r <= upper))
})
