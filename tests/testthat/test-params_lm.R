context("params_lm.R unit tests")
library("MASS")

# params_lm() works as expected ------------------------------------------------
test_that("params_lm works with matrix coefficients", {
  p <- params_lm(coefs = matrix(c(1, 2), nrow = 1),
                 sigma = 2)
  expect_equal(p$sigma, 2)
})

test_that("params_lm automatically names matrix columns", {
  p <- params_lm(coefs = matrix(c(1, 2), nrow = 1),
                 sigma = 2)
  expect_equal(colnames(p$coefs), c("x1", "x2"))
})

# params_lm() throws errors ----------------------------------------------------
test_that("params_lm throws error if numbers of samples are inconsistent ", {
  expect_error(
    params_lm(coefs = c(1, 3), sigma = 2),
    "Number of samples in 'sigma' is not equal to the number of samples in 'coefs'"
  )
})

test_that("params_lm throws error if sigma is not numeric", {
  expect_error(
    params_lm(coefs = matrix(c(1, 2), nrow = 1),
              sigma = "cat"),
    "is.numeric(sigma) is not TRUE",
    fixed = TRUE
  )
})

# create_params.lm() -----------------------------------------------------------
medcost_fit <- stats::lm(costs ~ female + state_name, 
                         data = psm4_exdata$costs$medical)

test_that("create_params.lm works with point estimates", {
  
  # point estimates
  p <- create_params(medcost_fit, n = 5, uncertainty = "none")
  expect_equal(p$coefs[, ], coef(medcost_fit))
  expect_equal(p$sigma, summary(medcost_fit)$sigma)
  expect_equal(colnames(p$coefs), names(medcost_fit$coefficients))
})

test_that("create_params.lm works with PSA", {
  set.seed(101)
  p <- create_params(medcost_fit, n = n, uncertainty = "normal")
  set.seed(101)
  r <- mvrnorm(n, medcost_fit$coefficients, vcov(medcost_fit))
  expect_equal(p$coefs, r)
})

test_that("create_params.lm works with PSA and sample size of 1", {
  p <- create_params(medcost_fit, n = 1, uncertainty = "normal")
  expect_error(pars_lm$coefs, NA)
})