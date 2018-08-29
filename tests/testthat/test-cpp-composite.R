context("integrate.cpp unit tests")
library("pracma")

f <- function(x) x * x
t <- seq(1, 15)
n.t <- length(t)
y <- f(t)

test_that("trapz",{
  expect_equal(hesim:::C_test_trapz(t, y), pracma::trapz(t, y))
})