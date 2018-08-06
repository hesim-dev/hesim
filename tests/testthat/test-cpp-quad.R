context("quad.h unit tests")

f <- function(x) {dnorm(x)}
test_that("definite_integral", {
  R_value <- stats::integrate(f, 0, 2)$value
  expect_equal(hesim:::test_quad_functor(0, 2), 
               R_value)

  expect_equal(hesim:::test_quad_lambda(0, 2),
               R_value)
})

test_that("indefinite_integral", {
  # Lower bound to infinity
  R_value <- stats::integrate(f, lower = 0, upper = Inf)$value
  expect_equal(hesim:::test_quad_functor(0, Inf),
               R_value)
  expect_equal(hesim:::test_quad_lambda(0, Inf),
               R_value)
  
  # -Infinity to upper bound
  R_value <- stats::integrate(f, lower = -Inf, upper = 1)$value
  expect_equal(hesim:::test_quad_lambda(-Inf, 1),
               R_value)
  
  # -Infinity to Infinity
  R_value <- stats::integrate(f, lower = -Inf, upper = Inf)$value
  expect_equal(hesim:::test_quad_lambda(-Inf, Inf),
               R_value)
})

