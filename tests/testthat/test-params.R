context("params.R unit tests")
library("flexsurv")

# Linear model -----------------------------------------------------------------
fit1 <- stats::lm(costs ~ female, data = psm4_exdata$costs$medical)
fit2 <- stats::lm(costs ~ 1, data = psm4_exdata$costs$medical)
fit_list <- lm_list(fit1 = fit1,
                    fit2 = fit2)

test_that("params_lm", {
  pars_lm <- params_lm(coefs = matrix(c(1, 2), nrow = 1),
                            sigma = 2)
  expect_equal(pars_lm$sigma, 2)
  
  # errors
  expect_error(pars_lm(coefs = c(1, 3), sigma = 2))
  expect_error(pars_lm(coefs = matrix(c(1, 2), nrow = 1),
                            sigma = "cat"))
  expect_error(pars_lm(coefs = matrix(c(1, 2), nrow = 1),
                            sigma = c(1, 2)))
})

test_that("create_params.lm", {
  n <- 5
  
  # errors
  expect_error(create_params(list(x = 5), 2))
  
  # point estimates
  pars_lm <- create_params(fit1, n = 5, uncertainty = "none")
  expect_equal(pars_lm$coefs[, ], coef(fit1))
  expect_equal(pars_lm$sigma, summary(fit1)$sigma)
  expect_equal(colnames(pars_lm$coefs), names(fit1$coefficients))
  
  # sampling
  set.seed(101)
  pars_lm <- create_params(fit1, n = n, uncertainty = "normal")
  set.seed(101)
  r <- MASS::mvrnorm(n, fit1$coefficients, vcov(fit1))
  expect_equal(pars_lm$coefs, r)
  
  # sample size of 1
  pars_lm <- create_params(fit1, n = 1, uncertainty = "normal")
  expect_error(pars_lm$coefs, NA)
})

test_that("params_lm_list", {
  pars_lm1 <- create_params(fit1, n = 5)
  pars_lm2 <- create_params(fit1, n = 5)
  pars_lm_list <- params_lm_list(pars_lm1, pars_lm2)
  expect_equal(length(pars_lm_list), 2)
  
  # errors
  pars_lm2 <- create_params(fit1, n = 2)
  expect_error(pars_lm_list(pars_lm1, pars_lm2))
})

test_that("create_params.lm_list", {
  pars_lm_list <- create_params(fit_list, n = 3)
  expect_equal(length(pars_lm_list), length(fit_list))
})
