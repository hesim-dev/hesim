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

test_that("form_params.lm", {
  n <- 5
  
  # errors
  expect_error(form_params(list(x = 5), 2))
  
  # point estimates
  pars_lm <- form_params(fit1, n = 5, point_estimate = TRUE)
  expect_equal(pars_lm$coefs[, ], coef(fit1))
  expect_equal(pars_lm$sigma, summary(fit1)$sigma)
  expect_equal(colnames(pars_lm$coefs), names(fit1$coefficients))
  
  # sampling
  set.seed(101)
  pars_lm <- form_params(fit1, n = n, point_estimate = FALSE)
  set.seed(101)
  r <- MASS::mvrnorm(n, fit1$coefficients, vcov(fit1))
  expect_equal(pars_lm$coefs, r)
  
  # sample size of 1
  pars_lm <- form_params(fit1, n = 1, point_estimate = FALSE)
  expect_error(pars_lm$coefs, NA)
})

test_that("params_lm_list", {
  pars_lm1 <- form_params(fit1, n = 5)
  pars_lm2 <- form_params(fit1, n = 5)
  pars_lm_list <- params_lm_list(pars_lm1, pars_lm2)
  expect_equal(length(pars_lm_list), 2)
  
  # errors
  pars_lm2 <- form_params(fit1, n = 2)
  expect_error(pars_lm_list(pars_lm1, pars_lm2))
})

test_that("form_params.lm_list", {
  pars_lm_list <- form_params(fit_list, n = 3)
  expect_equal(length(pars_lm_list), length(fit_list))
})

# Survival model ---------------------------------------------------------------
test_that("params_surv", {
  pars_surv <- params_surv(coefs = list(matrix(c(1, 2, 3, 4), nrow = 2)),
                            dist = "exponential")
  expect_equal(pars_surv$n_samples, 2)
  expect_true(inherits(pars_surv, "params_surv"))
  pars_surv <- params_surv(coefs = list(p1 = matrix(c(1, 2, 3, 4), nrow = 2),
                                          p2 = matrix(c(5, 6, 7, 8), nrow = 2)),
                            dist = "weibull")
  expect_equal(pars_surv$n_samples, 2)
  
  # errors
  expect_error(params_surv(coefs = matrix(c(1, 2, 3, 4), nrow = 2),
                            dist = "exponential"))
  expect_error(params_surv(coefs = list(rep(3, 3)),
                            dist = "exponential"))
  expect_error(params_surv(coefs = list(matrix(1), 3),
                            dist = "exponential"))
  expect_error(params_surv(coefs = list(matrix(c(1, 2), nrow = 1),
                                        matrix(c(1, 2, 3, 4), nrow = 2)),
                            dist = "weibull"))
})

test_that("form_params.flexsurv", {
  # no regressors
  ## exponential
  fit <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                     data = ovarian, dist = "exponential")
  pars_surv <- form_params(fit, point_estimate = TRUE)
  expect_equal(pars_surv$coefs$rate[, ], fit$res.t["rate", "est"])
  
  ### sample of size 1
  expect_error(form_params(fit, n = 1)$coefs$rate, NA)
  
  ## weibull
  fit <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                     data = ovarian, dist = "weibull")
  n <- 2
  set.seed(102)
  pars_surv <- form_params(fit, n = n)
  set.seed(102)
  sim <- flexsurv::normboot.flexsurvreg(fit, B = n, transform = TRUE)
  expect_equal(pars_surv$coefs$shape[, ], sim[, "shape"])
  expect_equal(pars_surv$coefs$scale[, ], sim[, "scale"])
  
  ## gengamma
  fit <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                     data = ovarian, dist = "gengamma")
  pars_surv <- form_params(fit)
  expect_equal(length(pars_surv$coefs), 3)
  
  ## covariates on 1 paramters
  fit <- flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                     data = ovarian, dist = "lognormal")
  pars_surv <- form_params(fit, n = 3)
  expect_equal(ncol(pars_surv$coefs$meanlog), 2)
  expect_equal(ncol(pars_surv$coefs$sdlog), 1)
  
  ## covariates on 2 paramters
  fit <- flexsurv::flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                     anc = list(sigma = ~ group), dist = "gengamma") 
  pars_surv <- form_params(fit, n = 2)
  expect_equal(ncol(pars_surv$coefs$mu), 3)
  expect_equal(ncol(pars_surv$coefs$sigma), 3)
  expect_equal(ncol(pars_surv$coefs$Q), 1)
  
  # spline
  fit <- flexsurv::flexsurvspline(Surv(recyrs, censrec) ~ group, data = bc, k = 1, 
                        scale = "hazard")
  pars_surv <- form_params(fit, n = 2)
  expect_error(pars_surv$coefs, NA)
})

fit_exp <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                                 data = ovarian, dist = "exp")
fit_wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                                 data = ovarian, dist = "weibull")
fit_lnorm <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                                 data = ovarian, dist = "lognormal")
pars_surv_exp <- form_params(fit_exp, n = 2)
pars_surv_wei <- form_params(fit_wei, n = 2)
pars_surv_lnorm <- form_params(fit_lnorm, n = 2)

fit_list1 <- flexsurvreg_list(exp = fit_exp, wei = fit_wei)
fit_list2 <- flexsurvreg_list(exp = fit_exp, wei = fit_lnorm)


test_that("params_surv_list", {
  pars_surv_list <- params_surv_list(pars_surv_exp, pars_surv_wei)
  expect_equal(length(pars_surv_list), 2)
})

test_that("form_params.flexsurvreg_list", {
  pars_surv_list <- form_params(fit_list1, n = 3)
  expect_true(inherits(pars_surv_list, "params_surv_list"))
})

test_that("form_params.partsurvfit", {
  fit1 <- flexsurvspline(Surv(endpoint1_time, endpoint1_status) ~ age,
                         data = psm4_exdata$survival[1:30, ])
  fit2 <- flexsurvspline(Surv(endpoint2_time, endpoint2_status) ~ age,
                         data = psm4_exdata$survival[1:30, ])
  partsurv_fit <- partsurvfit(flexsurvreg_list(fit1 = fit1, fit2 = fit2), data = psm4_exdata$survival)
  pars_surv_list <- form_params(partsurv_fit, n = 2)
  expect_equal(length(pars_surv_list), 2)
})

test_that("params_joined_surv", {
  pars_joined_surv <- params_joined_surv(exp = pars_surv_exp,
                                         wei = pars_surv_wei,
                                         times = 3)
  expect_true(inherits(pars_joined_surv, "params_joined_surv"))
  expect_equal(length(pars_joined_surv$models), 2)
  expect_equal(pars_joined_surv$times, 3)
  
  # errors
  expect_error(params_joined_surv(exp = pars_surv_exp, 5, times = 3))
})

test_that("params_joined_surv_list", {
  pars_surv_list1 <- params_surv_list(pars_surv_exp, pars_surv_wei)
  pars_surv_list2 <- params_surv_list(pars_surv_exp, pars_surv_lnorm)
  pars_joined_surv_list <- params_joined_surv_list(model1 = pars_surv_list1,
                                                     model2 = pars_surv_list2,
                                                     times = list(3, 5))
  expect_true(inherits(pars_joined_surv_list, "params_joined_surv_list"))
  expect_equal(length(pars_joined_surv_list$models), 2)
  
  # errors
  expect_error(params_joined_surv_list(pars_surv_exp, pars_surv_wei, times = 3))
  expect_error(params_joined_surv_list(model1 = pars_surv_list1,
                                       model2 = pars_surv_list2, 
                                       times = 3))
})

test_that("form_params.joined_flexsurvreg", {
  joined_fsreg <- joined_flexsurvreg(exp = fit_exp, wei = fit_wei,
                                           times = 6)
  pars_joined_surv <- form_params(joined_fsreg, n = 2)
  
  expect_true(inherits(pars_joined_surv, "params_joined_surv"))
  expect_equal(length(pars_joined_surv$models), 2)
})

test_that("form_params.joined_flexsurvreg_list", {
  joined_fsreg_list <- joined_flexsurvreg_list(fit1 = fit_list1, fit2 = fit_list2,
                                            times = list(1, 2))
  pars_joined_surv_list <- form_params(joined_fsreg_list, n = 2)
  
  expect_true(inherits(pars_joined_surv_list, "params_joined_surv_list"))
  expect_equal(length(pars_joined_surv_list$models), 2)
})