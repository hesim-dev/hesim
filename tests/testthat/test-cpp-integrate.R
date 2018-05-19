context("integrate.cpp unit tests")
library("pracma")

f <- function(x) x * x
t <- seq(1, 15)
n.t <- length(t)
y <- f(t)

test_that("trapz",{
  expect_equal(hesim:::C_test_trapzfun(t), pracma::trapz(t, y))
  expect_equal(hesim:::C_test_trapz(t, y), pracma::trapz(t, y))
  expect_equal(hesim:::C_test_trapzfun(t), 
               pracma::cotes(f, t[1], t[n.t], n = n.t - 1, nodes = 2))
})

test_that("cumtrapzfun",{
  expect_equal(hesim:::C_test_cumtrapzfun(t), c(pracma::cumtrapz(t, y)))
})

test_that("simpsfun",{
  expect_equal(hesim:::C_test_simpsfun(t),
               pracma::cotes(f, t[1], t[n.t], n = n.t - 1, nodes = 3))
})

test_that("cumsimpsfun",{
  pracma.cumsimps <- rep(NA, n.t)
  pracma.cumsimps[1] <- 0
  for (i in 1:(n.t - 1)){
    pracma.cumsimps[i + 1] <- pracma.cumsimps[i] +
        pracma::cotes(f, t[i], t[i + 1], n = n.t - 1, nodes = 3)
  }
  hesim.cumsimps <- hesim:::C_test_cumsimpsfun(t)
  expect_equal(pracma.cumsimps, hesim.cumsimps)
})