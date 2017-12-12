context("survival.R unit tests")
library("flexsurv")

# Sample parameters ------------------------------------------------------------
test_that("rsample.flexsurvreg", {
  
  # covariates on no paramters
  fit <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                     data = ovarian, dist = "weibull")
  n <- 5
  set.seed(101)
  fit.sample <- rsample.flexsurvreg(fit, n = n)
  set.seed(101)
  fit.normboot <- normboot.flexsurvreg(fit, B = n)
  expect_equal(length(fit.sample), length(fit$dlist$pars))
  expect_equal(c(fit.sample$shape), fit.normboot[, "shape"])
  expect_equal(nrow(fit.sample$scale), n)
  expect_false(any(lapply(fit.sample, ncol) != 1))
  
  ## sample of size 1
  expect_error(rsample.flexsurvreg(fit, n = 1), NA)
  
  ## check rsample()
  expect_error(rsample(fit, n = n), NA)
  
  # covariates on 1 paramters
  ## lognormal
  fit <- flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                       data = ovarian, dist = "lognormal")
  fit.sample <- rsample.flexsurvreg(fit, n = 3)
  expect_equal(ncol(fit.sample$meanlog), 2)
  expect_equal(ncol(fit.sample$sdlog), 1)
  
  # covariates on 2 parameters
  fit <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                       anc = list(sigma = ~ group), dist = "gengamma") 
  fit.sample <- rsample.flexsurvreg(fit, n = 2)
  expect_equal(ncol(fit.sample$mu), 3)
  expect_equal(ncol(fit.sample$sigma), 3)
  expect_equal(ncol(fit.sample$Q), 1)
})

test_that("sample.flexsurvspline", {
  fit <- flexsurvspline(Surv(recyrs, censrec) ~ group, data = bc, k=1, 
                        scale = "hazard")
  fit.normboot <- rsample(fit, n = 2)
  expect_error(fit.normboot, NA)


})

