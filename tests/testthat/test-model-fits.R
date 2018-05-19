context("model-fits unit tests")
library("flexsurv")

# formula ----------------------------------------------------------------------
test_that("joined_formula", {
  joined.formula <- joined_formula(f1 = formula(~ age), f2 = formula(~ 1),
                                   f3 = formula(~age + age2),
                                   times = c(2, 3))
  expect_true(inherits(joined.formula, "joined_formula"))
  
  # errors
  expect_error(joined_formula(f1 = formula(~ age), f2 = formula(~ 1),
                              f3 = formula(~age + age2),
                              times = c(3, 2)))
  expect_error(joined_formula(f1 = formula(~ age), f2 = formula(~ 1),
                              f3 = formula(~age + age2),
                              times = 3))
  expect_error(joined_formula(2, times = 3))
  expect_error(joined_formula(formula(~1), times = "times"))
})

test_that("formula_list", {
  f.list <- formula_list(f1 = formula(~ age), f2 = formula(~ 1))
  expect_true(inherits(f.list, "formula_list"))
  
  f.list <- formula_list(c(f1 = formula(~ age), f2 = formula(~ 1)))
  expect_true(inherits(f.list, "formula_list"))
  
  f.list <- formula_list(list(f1 = formula(~ age), f2 = formula(~ 1)))
  expect_true(inherits(f.list, "formula_list"))
  
  expect_error(formula_list(3))
  expect_error(formula_list(list(3)))
})

test_that("Nested formula_list", {
  f.list1 <- formula_list(f1 = formula(~ age), f2 = formula(~ 1))
  f.list2 <- formula_list(f1 = ~ female)
  f.list <- formula_list(f1 = f.list1, f2 = f.list2)
  expect_true(inherits(f.list, "formula_list"))
})

test_that("joined_formula_list", {
  f.list1 <- formula_list(f1 = formula(~ age), f2 = formula(~ 1))
  f.list2 <- formula_list(f1 = formula(~ age), f2 = formula(~ age2), f3 = formula(~age3))
  joined.formula.list <- joined_formula_list(model1 = f.list1, model2 = f.list2, 
                                             times = list(5, c(2, 4)))
  expect_true(inherits(joined.formula.list, "joined_formula_list"))
  expect_equal(joined.formula.list$times, list(5, c(2, 4)))
  joined.formula.list <- joined_formula_list(list(model1 = f.list1, model2 = f.list2), 
                                             times = list(5, c(2, 4)))
  expect_true(inherits(joined.formula.list, "joined_formula_list"))
  
  
  # errors
  expect_error(joined_formula_list(model1 = f.list1, model2 = f.list2,
                                   times = 5))
  expect_error(joined_formula_list(model1 = f.list1, model2 = f.list2,
                                   times = list(5, c(3, 2))))
})

# lm ---------------------------------------------------------------------------
dat <- part_surv4_simdata$costs$medical
fit1 <- stats::lm(costs ~ 1, data = dat)
fit2 <- stats::lm(costs ~ female, data = dat)
lm.list <- lm_list(fit1 = fit1, fit2 = fit2)

test_that("lm_list", {
  expect_true(inherits(lm.list, "lm_list"))

  lm.list <- lm_list(list(fit1 = fit1, fit2 = fit2))
  expect_true(inherits(lm.list, "lm_list"))
})

# flexsurvreg ------------------------------------------------------------------
dat <- part_surv4_simdata$survival
fit1a <- flexsurv::flexsurvreg(formula = Surv(endpoint1_time, endpoint1_status) ~ 1, 
                              data = dat, dist = "weibull")
fit1b <- flexsurv::flexsurvreg(formula = Surv(endpoint1_time, endpoint1_status) ~ 1, 
                              data = dat, dist = "exp")
flexsurvreg.list1 <- flexsurvreg_list(wei = fit1a, exp = fit1b)
  
fit2a <- flexsurv::flexsurvreg(formula = Surv(endpoint1_time, endpoint1_status) ~ age, 
                               data = dat, dist = "weibull")
fit2b <- flexsurv::flexsurvreg(formula = Surv(endpoint1_time, endpoint1_status) ~ age, 
                               data = dat, dist = "exp")
fit2c <- flexsurv::flexsurvreg(formula = Surv(endpoint1_time, endpoint1_status) ~ age, 
                               data = dat, dist = "gompertz")
flexsurvreg.list2 <-  flexsurvreg_list(wei = fit2a, exp = fit2b, gomp = fit2c)

test_that("flexsurvreg_list", {
  expect_true(inherits(flexsurvreg.list1, "flexsurvreg_list"))
})

test_that("joined_flexsurvreg", {
  joined.flexsurvreg <- joined_flexsurvreg(fit1a, fit1b, times = 5)
  expect_true(inherits(joined.flexsurvreg, "joined_flexsurvreg"))
})

test_that("joined_flexsurvreg_list", {
  joined.flexsurvreg.list <- joined_flexsurvreg_list(flexsurvreg.list1, flexsurvreg.list2,
                                                times = list(1, c(2, 5)))
  expect_true(inherits(joined.flexsurvreg.list, "joined_flexsurvreg_list"))
})

# partsurvfit ------------------------------------------------------------------
test_that("partsurvfit", {
  part.surv.fit <- partsurvfit(object = flexsurvreg.list1, data = ovarian)
  expect_true(inherits(part.surv.fit, "partsurvfit"))
})
