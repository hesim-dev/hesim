context("params.R unit tests")
library("flexsurv")

# Linear model -----------------------------------------------------------------
fit1 <- stats::lm(costs ~ female, data = part_surv4_simdata$costs$medical)
fit2 <- stats::lm(costs ~ 1, data = part_surv4_simdata$costs$medical)
fit.list <- lm_list(fit1 = fit1,
                    fit2 = fit2)

test_that("params_lm", {
  params.lm <- params_lm(coefs = matrix(c(1, 2), nrow = 1),
                            sigma = 2)
  expect_equal(params.lm$sigma, 2)
  
  # errors
  expect_error(params_lm(coefs = c(1, 3), sigma = 2))
  expect_error(params_lm(coefs = matrix(c(1, 2), nrow = 1),
                            sigma = "cat"))
  expect_error(params_lm(coefs = matrix(c(1, 2), nrow = 1),
                            sigma = c(1, 2)))
})

test_that("form_params.lm", {
  n <- 5
  
  # errors
  expect_error(form_params(list(x = 5), 2))
  
  # point estimates
  params.lm <- form_params(fit1, n = 5, point_estimate = TRUE)
  expect_equal(params.lm$coefs[, ], coef(fit1))
  expect_equal(params.lm$sigma, summary(fit1)$sigma)
  expect_equal(colnames(params.lm$coefs), names(fit1$coefficients))
  
  # sampling
  set.seed(101)
  params.lm <- form_params(fit1, n = n, point_estimate = FALSE)
  set.seed(101)
  r <- MASS::mvrnorm(n, fit1$coefficients, vcov(fit1))
  expect_equal(params.lm$coefs, r)
  
  # sample size of 1
  params.lm <- form_params(fit1, n = 1, point_estimate = FALSE)
  expect_error(params.lm$coefs, NA)
})

test_that("params_lm_list", {
  params.lm1 <- form_params(fit1, n = 5)
  params.lm2 <- form_params(fit1, n = 5)
  params.lm.list <- hesim:::params_lm_list(params.lm1, params.lm2)
  expect_equal(length(params.lm.list), 2)
  
  # errors
  params.lm2 <- form_params(fit1, n = 2)
  expect_error(hesim:::params_lm_list(params.lm1, params.lm2))
})

test_that("form_params.lm_list", {
  params.lm.list <- form_params(fit.list, n = 3)
  expect_equal(length(params.lm.list), length(fit.list))
})

# Survival model ---------------------------------------------------------------
test_that("params_surv", {
  params.surv <- params_surv(coefs = list(matrix(c(1, 2, 3, 4), nrow = 2)),
                            dist = "exponential")
  expect_equal(params.surv$n_samples, 2)
  expect_true(inherits(params.surv, "params_surv"))
  params.surv <- params_surv(coefs = list(p1 = matrix(c(1, 2, 3, 4), nrow = 2),
                                          p2 = matrix(c(5, 6, 7, 8), nrow = 2)),
                            dist = "weibull")
  expect_equal(params.surv$n_samples, 2)
  
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
  params.surv <- form_params(fit, point_estimate = TRUE)
  expect_equal(params.surv$coefs$rate[, ], fit$res.t["rate", "est"])
  
  ### sample of size 1
  expect_error(form_params(fit, n = 1)$coefs$rate, NA)
  
  ## weibull
  fit <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                     data = ovarian, dist = "weibull")
  n <- 2
  set.seed(102)
  params.surv <- form_params(fit, n = n)
  set.seed(102)
  sim <- flexsurv::normboot.flexsurvreg(fit, B = n, transform = TRUE)
  expect_equal(params.surv$coefs$shape[, ], sim[, "shape"])
  expect_equal(params.surv$coefs$scale[, ], sim[, "scale"])
  
  ## gengamma
  fit <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                     data = ovarian, dist = "gengamma")
  params.surv <- form_params(fit)
  expect_equal(length(params.surv$coefs), 3)
  
  ## covariates on 1 paramters
  fit <- flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                     data = ovarian, dist = "lognormal")
  params.surv <- form_params(fit, n = 3)
  expect_equal(ncol(params.surv$coefs$meanlog), 2)
  expect_equal(ncol(params.surv$coefs$sdlog), 1)
  
  ## covariates on 2 paramters
  fit <- flexsurv::flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                     anc = list(sigma = ~ group), dist = "gengamma") 
  params.surv <- form_params(fit, n = 2)
  expect_equal(ncol(params.surv$coefs$mu), 3)
  expect_equal(ncol(params.surv$coefs$sigma), 3)
  expect_equal(ncol(params.surv$coefs$Q), 1)
  
  # spline
  fit <- flexsurv::flexsurvspline(Surv(recyrs, censrec) ~ group, data = bc, k = 1, 
                        scale = "hazard")
  params.surv <- form_params(fit, n = 2)
  expect_error(params.surv$coefs, NA)
})

fit.exp <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                                 data = ovarian, dist = "exp")
fit.wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                                 data = ovarian, dist = "weibull")
fit.lnorm <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                                 data = ovarian, dist = "lognormal")
params.surv.exp <- form_params(fit.exp, n = 2)
params.surv.wei <- form_params(fit.wei, n = 2)
params.surv.lnorm <- form_params(fit.lnorm, n = 2)

fit.list1 <- flexsurvreg_list(exp = fit.exp, wei = fit.wei)
fit.list2 <- flexsurvreg_list(exp = fit.exp, wei = fit.lnorm)


test_that("params_surv_list", {
  params.surv.list <- params_surv_list(params.surv.exp, params.surv.wei)
  expect_equal(length(params.surv.list), 2)
})

test_that("form_params.flexsurvreg_list", {
  params.surv.list <- form_params(fit.list1, n = 3)
  expect_true(inherits(params.surv.list, "params_surv_list"))
})

test_that("form_params.partsurvfit", {
  fit1 <- flexsurvspline(Surv(endpoint1_time, endpoint1_status) ~ age,
                         data = part_surv4_simdata$survival[1:30, ])
  fit2 <- flexsurvspline(Surv(endpoint2_time, endpoint2_status) ~ age,
                         data = part_surv4_simdata$survival[1:30, ])
  partsurv.fit <- partsurvfit(flexsurvreg_list(fit1 = fit1, fit2 = fit2), data = part_surv4_simdata$survival)
  params.surv.list <- form_params(partsurv.fit, n = 2)
  expect_equal(length(params.surv.list), 2)
})

test_that("params_joined_surv", {
  params.joined.surv <- params_joined_surv(exp = params.surv.exp,
                                           wei = params.surv.wei,
                                           times = 3)
  expect_true(inherits(params.joined.surv, "params_joined_surv"))
  expect_equal(length(params.joined.surv$models), 2)
  expect_equal(params.joined.surv$times, 3)
  
  # errors
  expect_error(params_joined_surv(exp = params.surv.exp, 5, times = 3))
})

test_that("params_joined_surv_list", {
  params.surv.list1 <- params_surv_list(params.surv.exp, params.surv.wei)
  params.surv.list2 <- params_surv_list(params.surv.exp, params.surv.lnorm)
  params.joined.surv.list <- params_joined_surv_list(model1 = params.surv.list1,
                                                     model2 = params.surv.list2,
                                                     times = list(3, 5))
  expect_true(inherits(params.joined.surv.list, "params_joined_surv_list"))
  expect_equal(length(params.joined.surv.list$models), 2)
  
  # errors
  expect_error(params_joined_surv_list(params.surv.exp, params.surv.wei, times = 3))
  expect_error(params_joined_surv_list(model1 = params.surv.list1,
                                       model2 = params.surv.list2, 
                                       times = 3))
})

test_that("form_params.joined_flexsurvreg", {
  joined.flexsurvreg <- joined_flexsurvreg(exp = fit.exp, wei = fit.wei,
                                           times = 6)
  params.joined.surv <- form_params(joined.flexsurvreg, n = 2)
  
  expect_true(inherits(params.joined.surv, "params_joined_surv"))
  expect_equal(length(params.joined.surv$models), 2)
})

test_that("form_params.joined_flexsurvreg_list", {
  joined.flexsurvreg.list <- joined_flexsurvreg_list(fit1 = fit.list1, fit2 = fit.list2,
                                                times = list(1, 2))
  params.joined.surv.list <- form_params(joined.flexsurvreg.list, n = 2)
  
  expect_true(inherits(params.joined.surv.list, "params_joined_surv_list"))
  expect_equal(length(params.joined.surv.list$models), 2)
})