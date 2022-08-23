context("params_lm.R unit tests")
library("MASS")

p_ex <- params_lm(
  coefs = matrix(c(1, 2), nrow = 1),
  sigma = 2
)

# params_lm() works as expected ------------------------------------------------
test_that("params_lm works with matrix coefficients", {
  expect_equal(p_ex$sigma, 2)
})

test_that("params_lm automatically names matrix columns", {
  expect_equal(colnames(p_ex$coefs), c("x1", "x2"))
})

test_that("params_lm works with data.frame coefficients", {
  p <- params_lm(
    coefs = data.frame(intercept = c(1, 2)),
    sigma = 2
  )
  expect_equal(p$coefs[, "intercept"], c(1, 2))
})

test_that("params_lm works with default sigma", {
  p <- params_lm(coefs = data.frame(intercept = c(1, 2)))
  expect_equal(p$sigma, c(1, 1))
})

# params_lm() throws errors ----------------------------------------------------
test_that("params_lm throws error if numbers of samples are inconsistent ", {
  expect_error(
    params_lm(coefs = c(1, 3, 2), sigma = c(2, 2)),
    "Number of samples in 'sigma' is not equal to the number of samples in 'coefs'"
  )
})

test_that("params_lm throws error if sigma is not numeric", {
  expect_error(
    params_lm(
      coefs = matrix(c(1, 2), nrow = 1),
      sigma = "cat"
    ),
    "is.numeric(sigma) is not TRUE",
    fixed = TRUE
  )
})

# summary_params.lm() ----------------------------------------------------------
test_that("summary.params_lm()", {
  ps <- summary(p_ex)
  expect_equal(ps$term, c("x1", "x2", "sigma"))
  expect_equal(ps$parameter, c("mean", "mean", "sd"))
  expect_equal(
    unname(unlist(ps[term == "x1", .(mean, `2.5%`, `97.5%`)])),
    rep(1, 3)
  )
})

# print_params.lm() ------------------------------------------------------------
test_that("print.params_lm()", {
  expect_output(print(p_ex), "A \"params_lm\" object")
  expect_output(print(p_ex), "Summary of coefficients:")
  expect_output(print(p_ex), "Summary of sigma:")
})

# create_params.lm() -----------------------------------------------------------
medcost_fit <- lm(costs ~ female + state_name,
  data = psm4_exdata$costs$medical
)

test_that("create_params.lm works with point estimates", {

  # point estimates
  p <- create_params(medcost_fit, n = 5, uncertainty = "none")
  expect_equal(p$coefs[, ], coef(medcost_fit))
  expect_equal(p$sigma, summary(medcost_fit)$sigma)
  expect_equal(colnames(p$coefs), names(medcost_fit$coefficients))
})

test_that("create_params.lm works with PSA", {
  set.seed(101)
  p <- create_params(medcost_fit, n = 2, uncertainty = "normal")
  set.seed(101)
  r <- mvrnorm(n = 2, medcost_fit$coefficients, vcov(medcost_fit))
  expect_equal(p$coefs, r)
})

test_that("create_params.lm works with PSA and sample size of 1", {
  p <- create_params(medcost_fit, n = 1, uncertainty = "normal")
  expect_true(inherits(p, "params_lm"))
})
