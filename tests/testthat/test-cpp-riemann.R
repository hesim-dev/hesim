context("riemann.h unit tests")

cum_riemann <- function(x, f){
  sum <- rep(NA, length(x))
  sum[1] <- 0
  for (i in 2:length(x)){
    step <- x[i] - x[i - 1]
    mid <- x[i - 1] + step/2
    sum[i] <- step * f(mid) + sum[i - 1]
  }
  return(sum)
}
f <- function(x) x^2

test_that("Riemann sum", {
  x <- seq(0, 1, .01)
  R_value <- cum_riemann(x, f)[length(x)]
  stats::integrate(f, 0, 1)
  expect_equal(hesim:::test_riemann_x2(x), 
               R_value)
})

test_that("Riemann sum", {
  x <- seq(0, 1, .01)
  R_value <- cum_riemann(x, f)
  stats::integrate(f, 0, 1)
  expect_equal(hesim:::test_cum_riemann_x2(x), 
               R_value)
})

