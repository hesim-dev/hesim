context("model-fits unit tests")
library("flexsurv")

# formula ----------------------------------------------------------------------
test_that("formula_list", {
  f_list <- formula_list(f1 = formula(~age), f2 = formula(~1))
  expect_true(inherits(f_list, "formula_list"))

  f.list <- formula_list(c(f1 = formula(~age), f2 = formula(~1)))
  expect_true(inherits(f_list, "formula_list"))

  f.list <- formula_list(list(f1 = formula(~age), f2 = formula(~1)))
  expect_true(inherits(f_list, "formula_list"))

  expect_error(formula_list(3))
  expect_error(formula_list(list(3)))
})

test_that("Nested formula_list", {
  f_list1 <- formula_list(f1 = formula(~age), f2 = formula(~1))
  f_list2 <- formula_list(f1 = ~female)
  f_list <- formula_list(f1 = f_list1, f2 = f_list2)
  expect_true(inherits(f_list, "formula_list"))
})

# lm ---------------------------------------------------------------------------
dat <- psm4_exdata$costs$medical
fit1 <- stats::lm(costs ~ 1, data = dat)
fit2 <- stats::lm(costs ~ female, data = dat)
lm_list <- lm_list(fit1 = fit1, fit2 = fit2)

test_that("lm_list", {
  expect_true(inherits(lm_list, "lm_list"))

  lm_list <- lm_list(list(fit1 = fit1, fit2 = fit2))
  expect_true(inherits(lm_list, "lm_list"))
})

# flexsurvreg ------------------------------------------------------------------
dat <- psm4_exdata$survival
fit1a <- flexsurv::flexsurvreg(
  formula = Surv(endpoint1_time, endpoint1_status) ~ 1,
  data = dat, dist = "weibull"
)
fit1b <- flexsurv::flexsurvreg(
  formula = Surv(endpoint1_time, endpoint1_status) ~ 1,
  data = dat, dist = "exp"
)
flexsurvreg_list1 <- flexsurvreg_list(wei = fit1a, exp = fit1b)

fit2a <- flexsurv::flexsurvreg(
  formula = Surv(endpoint1_time, endpoint1_status) ~ age,
  data = dat, dist = "weibull"
)
fit2b <- flexsurv::flexsurvreg(
  formula = Surv(endpoint1_time, endpoint1_status) ~ age,
  data = dat, dist = "exp"
)
fit2c <- flexsurv::flexsurvreg(
  formula = Surv(endpoint1_time, endpoint1_status) ~ age,
  data = dat, dist = "gompertz"
)
flexsurvreg_list2 <- flexsurvreg_list(wei = fit2a, exp = fit2b, gomp = fit2c)

test_that("flexsurvreg_list", {
  expect_true(inherits(flexsurvreg_list1, "flexsurvreg_list"))
})

# partsurvfit ------------------------------------------------------------------
test_that("partsurvfit", {
  fit <- partsurvfit(object = flexsurvreg_list1, data = ovarian)
  expect_true(inherits(fit, "partsurvfit"))
})
