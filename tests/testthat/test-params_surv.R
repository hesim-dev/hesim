context("params_surv.R unit tests")
library("flexsurv")

# params_surv() works as expected ----------------------------------------------
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

# params_surv() throws errors --------------------------------------------------
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

test_that("params_surv() throws error if knots are not specified for a spline model", {
  expect_error(
    params_surv(coefs = list(matrix(.5)),
                aux = list(scale = "log_cumhazard"),
                dist = "survspline"),
    "'knots' must be specified in a spline model."
  )
})

test_that("params_surv() throws error if hazard scale is wong for spline model", {
  choices <- c("log_cumhazard", "log_hazard", "log_cumodds", "inv_normal")
  expect_error(
    params_surv(coefs = list(gamma0 = matrix(.5),
                             gamma1 = matrix(0)),
                aux = list(knots = c(0, 10),
                           scale = "log"),
                dist = "survspline"),
    paste0("The auxiliary argument 'scale' must be one of ", 
           paste(dQuote(choices), collapse = ", "))
  )
})

test_that("params_surv() throws error if time scale is wong for spline model", {
  expect_error(
    params_surv(coefs = list(gamma0 = matrix(.5),
                             gamma1 = matrix(0)),
                aux = list(knots = c(0, 10),
                           scale = "log_hazard",
                           timescale = "wrong"),
                dist = "survspline"),
    paste0("The auxiliary argument 'timescale' must be one of ", 
           paste(dQuote(c("log", "identity")), collapse = ", "))
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

test_that("params_surv() throws error if numbers of parameters in fractional polynomial model is wrong", {
  expect_error(
    params_surv(coefs = list(matrix(.8),
                             matrix(.9)),
                aux = list(powers = c(-2, -1)),
                dist = "fracpoly"),
    paste0("The number of parameters in a fractional polynomial model must equal ", 
           "the number of powers plus 1.")
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

# print.params_surv() ----------------------------------------------------------
test_that("print.params_surv() works as expected", {
  p <- params_surv(coefs = list(rate = rep(3, 10)),
                   dist = "exp")
  expect_output(print(p), "A \"params_surv\" object")
  expect_output(print(p), "Summary of coefficient estimates:")
  expect_output(print(p), "Number of parameter samples: 10")
  expect_output(print(p), "Distribution: exp")
})

test_that("print.params_surv() works with piecewise exponential model", {
  p <- params_surv(coefs = list(rate1 = 1, rate = 2),
                   dist = "pwexp",
                   aux = list(times = c(1, 5)))
  expect_output(print(p), "Times: 1 5")
})

test_that("print.params_surv() works with survival splines", {
  p <- params_surv(coefs = list(gamma0 = 1, gamma1 = 2),
                   dist = "survspline",
                   aux = list(knots = c(1, 3)))
  expect_output(print(p), "Knots: 1 3")
  expect_output(print(p), "Scale: log_cumhazard")
  expect_output(print(p), "Time scale: log")
})

test_that("print.params_surv() works with fractional polynomials", {
  p <- params_surv(coefs = list(gamma0 = 1, gamma2 = 2),
                   dist = "fracpoly",
                   aux = list(powers = 1))
  expect_output(print(p), "Distribution: fracpoly")
  expect_output(print(p), "Powers: 1")
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