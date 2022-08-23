context("quad.h unit tests")

f <- function(x) dnorm(x)
test_that("Definite integral", {
  R_value <- stats::integrate(f, 0, 2)$value
  expect_equal(
    hesim:::test_quad_dnorm(0, 2),
    R_value
  )
})

test_that("Indefinite integral", {
  # Lower bound to infinity
  R_value <- stats::integrate(f, lower = 0, upper = Inf)$value
  expect_equal(
    hesim:::test_quad_dnorm(0, Inf),
    R_value
  )

  # -Infinity to upper bound
  R_value <- stats::integrate(f, lower = -Inf, upper = 1)$value
  expect_equal(
    hesim:::test_quad_dnorm(-Inf, 1),
    R_value
  )

  # -Infinity to Infinity
  R_value <- stats::integrate(f, lower = -Inf, upper = Inf)$value
  expect_equal(
    hesim:::test_quad_dnorm(-Inf, Inf),
    R_value
  )
})

test_that("Warnings and errors", {
  expect_warning(hesim:::test_quad_ier1())
  expect_warning(hesim:::test_quad_ier4())
  expect_warning(hesim:::test_quad_ier5())
})
