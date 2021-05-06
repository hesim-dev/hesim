context("params_surv_list.R unit tests")
library("flexsurv")

# params_surv_list() works as expected -----------------------------------------
test_that("params_surv_list() works as expected", {
  n <- 5
  p <- params_surv_list(
    # Model 1
    mod1 = params_surv(
      coefs = list(rate = data.frame(intercept = rlnorm(n, 0, 1))),
      dist = "exp"
    ),
    
    # Model 2
    mod2 = params_surv(
      coefs = list(rate = data.frame(intercept = rlnorm(n, 1, 1))),
      dist = "exp"
    )
  )
  
  expect_true(inherits(p, "params_surv_list"))
  expect_equal(length(p), 2)
  expect_equal(p[[1]]$n_samples, n)
})

# params_surv_list() throws errors ---------------------------------------------
test_that("params_surv_list() throws errors if the number of samples are inconsistent", {
  expect_error(
    params_surv_list(
      params_surv(
        coefs = list(rate = data.frame(intercept = c(2, 3))),
        dist = "exp"
      ),
      params_surv(
        coefs = list(rate = data.frame(intercept = 1)),
        dist = "exp"
      )
    ),
    "The number of samples in each 'params_surv' object must be the same."
  )
})

# create_params.flexsurvreg_list() ---------------------------------------------
fit_exp <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                                 data = ovarian, dist = "exp")
fit_wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                                 data = ovarian, dist = "weibull")
fit_lnorm <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                                   data = ovarian, dist = "lognormal")
params_surv_exp <- create_params(fit_exp, n = 2)
params_surv_wei <- create_params(fit_wei, n = 2)
params_surv_lnorm <- create_params(fit_lnorm, n = 2)

fit_list1 <- flexsurvreg_list(exp = fit_exp, wei = fit_wei)
fit_list2 <- flexsurvreg_list(exp = fit_exp, wei = fit_lnorm)


test_that("create_params.flexsurvreg_list() works as expected", {
  p <- create_params(fit_list1, n = 3)
  expect_true(inherits(p, "params_surv_list"))
})

# create_params.partsurvfit() --------------------------------------------------
test_that("create_params.partsurvfit() works as expected", {
  fit1 <- flexsurvspline(Surv(endpoint1_time, endpoint1_status) ~ age,
                         data = psm4_exdata$survival[1:30, ])
  fit2 <- flexsurvspline(Surv(endpoint2_time, endpoint2_status) ~ age,
                         data = psm4_exdata$survival[1:30, ])
  partsurv_fit <- partsurvfit(flexsurvreg_list(fit1 = fit1, fit2 = fit2), 
                              data = psm4_exdata$survival)
  pars_surv_list <- create_params(partsurv_fit, n = 2)
  expect_equal(length(pars_surv_list), 2)
})

# params_joined_surv() ---------------------------------------------------------
test_that("params_joined_surv", {
  pars_joined_surv <- params_joined_surv(exp = params_surv_exp,
                                         wei = params_surv_wei,
                                         times = 3)
  expect_true(inherits(pars_joined_surv, "params_joined_surv"))
  expect_equal(length(pars_joined_surv$models), 2)
  expect_equal(pars_joined_surv$times, 3)
  
  # errors
  expect_error(params_joined_surv(exp = params_surv_exp, 5, times = 3))
})