context("params_)surv.R unit tests")
library("flexsurv")

# params_surv() ----------------------------------------------------------------
test_that("params_surv() works as expected for various distributions", {
  ## exponential
  p <- params_surv(coefs = list(matrix(c(1, 2, 3, 4), nrow = 2)),
                  dist = "exponential")
  expect_equal(p$n_samples, 2)
  expect_true(inherits(p, "params_surv"))
  
  ## weibull
  p <- params_surv(coefs = list(p1 = matrix(c(1, 2, 3, 4), nrow = 2),
                                p2 = matrix(c(5, 6, 7, 8), nrow = 2)),
                  dist = "weibull")
  expect_equal(p$n_samples, 2)
})

test_that("params_surv() with auxillary arguments", {
  p <- params_surv(coefs = list(matrix(.8),
                                matrix(.9)),
                   aux = list(time = c(1, 2)),
                   dist = "pwexp")
  expect_equal(p$dist, "pwexp")
  expect_equal(p$aux$time, c(1, 2))
})

test_that("params_surv() with data.frame passed to coefs", {
  p <- params_surv(coefs = list(rate = data.frame(intercept = 1)),
                   dist = "exp")
  expect_true(inherits(p, "params_surv"))
  
  p <- params_surv(
    coefs = list(
      shape = data.frame(
        intercept = c(1, 2)),
      scale = data.frame(
        intercept = c(1, 3),
        var = c(1, 1))
      ),
    dist = "weibull"
  )
  expect_equal(ncol(p$coefs$scale), 2)
})

test_that("params_surv() with vector passed to coefs", {
  p <- params_surv(coefs = list(rate = rep(3, 10)),
                   dist = "exp")
  expect_equal(nrow(p$coefs$rate), 10)
})

test_that("params_surv() throws error if coef argument is not a list", {
  expect_error(
    params_surv(coefs = matrix(c(1, 2, 3, 4), nrow = 2),
                           dist = "exponential"),
    "'coefs' must be a list."
  )
})

test_that("params_surv() throws error if number of rows in coef matrices are unequal", {
  expect_error(
    params_surv(coefs = list(matrix(c(1, 2), nrow = 1),
                             matrix(c(1, 2, 3, 4), nrow = 2)),
                dist = "weibull"),
    "Number of rows in all 'coefs' matrices must be equal."
    )
})

test_that("params_surv() throws error if piecewise exponential if times aren't consistent with rates", {
  expect_error(
    params_surv(coefs = list(matrix(.8),
                            matrix(.9)),
                aux = list(time = c(1)),
                dist = "pwexp"),
    "The length of 'time' must equal the length of 'coefs'."
  )
})

# summary.params_surv() --------------------------------------------------------
test_that("summary.params_surv()", {
  p <- params_surv(
    coefs = list(
      shape = data.frame(
        intercept = c(1, 2)),
      scale = data.frame(
        intercept = c(1, 3),
        var = c(1, 1))
    ),
    dist = "weibull"
  )
  
  ps <- summary(p)
  expect_true(inherits(ps, "data.table"))
  expect_equal(ps$parameter, c("shape", "scale", "scale"))
  expect_equal(ps$term, c("intercept", "intercept", "var"))
  expect_equal(ps$estimate, c(1.5, 2, 1))
})

# create_params.flexsurv() -----------------------------------------------------
test_that("create_params.flexsurv()", {
  # no regressors
  ## exponential
  fit <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                               data = ovarian, dist = "exponential")
  pars_surv <- create_params(fit, uncertainty = "none")
  expect_equal(pars_surv$coefs$rate[, ], fit$res.t["rate", "est"])
  
  ### sample of size 1
  expect_error(create_params(fit, n = 1)$coefs$rate, NA)
  
  ## weibull
  fit <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                               data = ovarian, dist = "weibull")
  n <- 2
  set.seed(102)
  pars_surv <- create_params(fit, n = n)
  set.seed(102)
  sim <- flexsurv::normboot.flexsurvreg(fit, B = n, transform = TRUE)
  expect_equal(pars_surv$coefs$shape[, ], sim[, "shape"])
  expect_equal(pars_surv$coefs$scale[, ], sim[, "scale"])
  
  ## gengamma
  fit <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                               data = ovarian, dist = "gengamma")
  pars_surv <- create_params(fit)
  expect_equal(length(pars_surv$coefs), 3)
  
  ## covariates on 1 paramters
  fit <- flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                     data = ovarian, dist = "lognormal")
  pars_surv <- create_params(fit, n = 3)
  expect_equal(ncol(pars_surv$coefs$meanlog), 2)
  expect_equal(ncol(pars_surv$coefs$sdlog), 1)
  
  ## covariates on 2 paramters
  fit <- flexsurv::flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                               anc = list(sigma = ~ group), dist = "gengamma") 
  pars_surv <- create_params(fit, n = 2)
  expect_equal(ncol(pars_surv$coefs$mu), 3)
  expect_equal(ncol(pars_surv$coefs$sigma), 3)
  expect_equal(ncol(pars_surv$coefs$Q), 1)
  
  # spline
  fit <- flexsurv::flexsurvspline(Surv(recyrs, censrec) ~ group, data = bc, k = 1, 
                                  scale = "hazard")
  pars_surv <- create_params(fit, n = 2)
  expect_error(pars_surv$coefs, NA)
})

fit_exp <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                                 data = ovarian, dist = "exp")
fit_wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                                 data = ovarian, dist = "weibull")
fit_lnorm <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                                   data = ovarian, dist = "lognormal")
pars_surv_exp <- create_params(fit_exp, n = 2)
pars_surv_wei <- create_params(fit_wei, n = 2)
pars_surv_lnorm <- create_params(fit_lnorm, n = 2)

fit_list1 <- flexsurvreg_list(exp = fit_exp, wei = fit_wei)
fit_list2 <- flexsurvreg_list(exp = fit_exp, wei = fit_lnorm)

# params_surv_list() -----------------------------------------------------------
test_that("params_surv_list", {
  pars_surv_list <- params_surv_list(pars_surv_exp, pars_surv_wei)
  expect_equal(length(pars_surv_list), 2)
})

# create_params.flexsurvreg_list() ---------------------------------------------
test_that("create_params.flexsurvreg_list", {
  pars_surv_list <- create_params(fit_list1, n = 3)
  expect_true(inherits(pars_surv_list, "params_surv_list"))
})

# create_params.partsurvfit() --------------------------------------------------
test_that("create_params.partsurvfit", {
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
  pars_joined_surv <- params_joined_surv(exp = pars_surv_exp,
                                         wei = pars_surv_wei,
                                         times = 3)
  expect_true(inherits(pars_joined_surv, "params_joined_surv"))
  expect_equal(length(pars_joined_surv$models), 2)
  expect_equal(pars_joined_surv$times, 3)
  
  # errors
  expect_error(params_joined_surv(exp = pars_surv_exp, 5, times = 3))
})

# params_joined_surv_list() ----------------------------------------------------
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

# create_params.joined_flexsurvreg() -------------------------------------------
test_that("create_params.joined_flexsurvreg", {
  joined_fsreg <- joined_flexsurvreg(exp = fit_exp, wei = fit_wei,
                                     times = 6)
  pars_joined_surv <- create_params(joined_fsreg, n = 2)
  
  expect_true(inherits(pars_joined_surv, "params_joined_surv"))
  expect_equal(length(pars_joined_surv$models), 2)
})

# create_params.joined_flexsurvreg_list() --------------------------------------
test_that("create_params.joined_flexsurvreg_list", {
  joined_fsreg_list <- joined_flexsurvreg_list(fit1 = fit_list1, fit2 = fit_list2,
                                               times = list(1, 2))
  pars_joined_surv_list <- create_params(joined_fsreg_list, n = 2)
  
  expect_true(inherits(pars_joined_surv_list, "params_joined_surv_list"))
  expect_equal(length(pars_joined_surv_list$models), 2)
})