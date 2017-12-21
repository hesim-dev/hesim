context("survival.R unit tests")
library("flexsurv")

# Sample parameters ------------------------------------------------------------
test_that("rsample.flexsurvreg", {
  
  # error when not an object of class flexsurvreg 
  expect_error(rsample(list(x = 5), 2))
  expect_error(rsample(matrix(1), n = 2))
  
  # covariates on no parameters
  fit <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                     data = ovarian, dist = "weibull")
  n <- 5
  set.seed(101)
  fit.sample <- rsample(fit, n = n)
  set.seed(101)
  fit.normboot <- normboot.flexsurvreg(fit, B = n)
  expect_equal(length(fit.sample$coefs), length(fit$dlist$pars))
  expect_equal(c(fit.sample$coefs$shape), fit.normboot[, "shape"])
  expect_equal(nrow(fit.sample$coefs$scale), n)
  expect_false(any(lapply(fit.sample$coefs[1:2], ncol) != 1))
  expect_true(inherits(fit.sample, "surv_pars"))
  
  ## sample of size 1
  expect_error(rsample(fit, n = 1), NA)
  
  # covariates on 1 paramters
  ## lognormal
  fit <- flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                     data = ovarian, dist = "lognormal")
  fit.sample <- rsample(fit, n = 3)
  expect_equal(ncol(fit.sample$coefs$meanlog), 2)
  expect_equal(ncol(fit.sample$coefs$sdlog), 1)
  
  # covariates on 2 parameters
  fit <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                     anc = list(sigma = ~ group), dist = "gengamma") 
  fit.sample <- rsample(fit, n = 2)
  expect_equal(ncol(fit.sample$coefs$mu), 3)
  expect_equal(ncol(fit.sample$coefs$sigma), 3)
  expect_equal(ncol(fit.sample$coefs$Q), 1)
  
  # spline
  fit <- flexsurvspline(Surv(recyrs, censrec) ~ group, data = bc, k=1, 
                        scale = "hazard")
  fit.sample <- rsample(fit, n = 2)
  expect_error(fit.sample, NA)
  
  # multiple models
  fit1 <- flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                      data = ovarian, dist = "exponential")
  fit2 <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                      data = ovarian, dist = "gamma")
  fits <- list(f1 = fit1, f2 = fit2)
  fits.sample <- rsample(fits, n = 2)
  expect_equal(fits.sample[[1]]$dist, "exp")
  expect_equal(fits.sample[[2]]$dist, "gamma")
  expect_true(inherits(fits.sample, "survlist_pars"))
  expect_error(rsample(list(x = 1), n = 1))
})