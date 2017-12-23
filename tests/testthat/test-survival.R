context("survival.R unit tests")
library("flexsurv")

# Sample parameters ------------------------------------------------------------
test_that("rsample", {
  
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
  fit.normboot <- normboot.flexsurvreg(fit, B = n, transform = TRUE)
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

# Extract design matrices ------------------------------------------------------
test_that("design_matrix", {
  # exponential - no covariates
  fit1 <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                     data = ovarian, dist = "exponential")
  dm.rate <- design_matrix(fit1)$rate
  expect_true(all.equal(c(dm.rate), rep(1, nrow(dm.rate))))
  
  # gompertz - covariate on rate parameter
  fit2 <- flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                     data = ovarian, dist = "gompertz")
  dm <- design_matrix(fit2)
  expect_equal(dm$rate[1, "age"], fit2$data$m$age[1])
  expect_equal(ncol(dm$rate), 2)
  expect_equal(length(dm), 2)
  fit2.samp <- rsample(fit2, n = 2)
  expect_equal(names(fit2.samp$coefs), names(dm))
  
  # list
  fit.list <- list(fit1, fit2)
  dm <- design_matrix(fit.list)
  expect_equal(length(dm), 2)
  expect_true(names(dm[[1]]) == "rate")
  expect_equal(dm[[2]]$rate[1, "age"], fit2$data$m$age[1])
})

# Survival R6 class ------------------------------------------------------------
test_that("Survival class no regressors", {
  fit <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                     data = ovarian, dist = "exponential")
  surv <- Survival$new(dist_name = fit$dlist$name,
                       surv_coefs = rsample(fit, n = 5)$coefs,
                       surv_X = list(fit$data$mml$rate))
  
  # Quantiles
  surv.quantiles <- surv$quantiles(q = c(.5, .8))
  expect_equal(stats::qexp(.5, exp(surv$surv_coefs$rate[1])), 
               surv.quantiles[sim == 1 & id == 1 & q == .5, value])
  
  # Survival curves
  surv.sc <- surv$surv(t = c(0, 1, 2))
  expect_equal(1 - stats::pexp(1, rate = exp(surv$surv_coefs$rate[2])),
                surv.sc[sim == 2 & id == 1 & t == 1, value])
  
  # Cummulative hazards
  surv.cumhaz <- surv$cumhazard(t = c(0, 1, 2))
  expect_equal(flexsurv::Hexp(2, rate = exp(surv$surv_coefs$rate[3])),
               surv.cumhaz[sim == 3 & id == 1 & t == 2, value])
  
  # Hazards
  surv.haz <- surv$hazard(t = c(0, 1, 2, 3))
  expect_equal(flexsurv::hexp(3, rate = exp(surv$surv_coefs$rate[3])),
               surv.haz[sim == 3 & id == 1 & t == 3, value])
  
  # rmst
  expect_true(surv$rmst(t = c(15, 20), r = 0)[1, "value"] > 
                surv$rmst(t = c(15, 20), r = 0.03)[1, "value"])
})

test_that("Survival class with regressors", {
  fit <- flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                     data = ovarian, dist = "weibull")
  surv <- Survival$new(dist_name = fit$dlist$name,
                       surv_coefs = rsample(fit, n = 5)$coefs,
                       surv_X = list(shape = matrix(1, nrow = nrow(fit$data$mml$scale)),
                                     scale = fit$data$mml$scale))
  
  # Quantiles
  surv.quantiles <- surv$quantiles(q = c(.5))
  scale1 <- surv$surv_X$scale[1,, drop = FALSE] %*% t(surv$surv_coefs$scale[1,, drop = FALSE])
  expect_equal(stats::qweibull(.5, shape = exp(surv$surv_coefs$shape[1]), 
                               scale = exp(scale1)), 
               surv.quantiles[sim == 1 & id == 1 & q == .5, value])
  
  # Survival curves
  surv.sc <- surv$surv(t = c(100, 365))
  expect_equal(1 - stats::pweibull(365, shape = exp(surv$surv_coefs$shape[1]), 
                               scale = exp(scale1)), 
               surv.sc[sim == 1 & id == 1 & t == 365, value])
})

test_that("Compare survival class with flexsurv", {
  flexsurv_comp <- function(dist_name){
      fit <- flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                        data = ovarian, dist = dist_name)
      surv <- Survival$new(dist_name = fit$dlist$name,
                           surv_coefs = rsample(fit, point_estimate = TRUE)$coefs,
                           surv_X = design_matrix(fit, indices = 1))
      
      # survival probabilities
      flexsurv.sc <- summary(fit, newdata = data.frame(age = fit$data$m$age[1]),
                             t = 100)
      surv.sc <- surv$surv(t = 100)
      expect_equal(flexsurv.sc[[1]]$est, surv.sc$value)
      
      # restricted mean survival time
      flexsurv.rmst <- summary(fit, newdata = data.frame(age = fit$data$m$age[1]),
                             t = 90, type = "rmst")
      surv.rmst <- surv$rmst(t = 90, r = 0)
      expect_equal(flexsurv.rmst[[1]]$est, surv.rmst$value, scale = 1, 
                   tol = .01)
  }
  flexsurv_comp("exp")
  flexsurv_comp("exponential")
  flexsurv_comp("weibull")
  flexsurv_comp("gamma")
  flexsurv_comp("gompertz")
  flexsurv_comp("llogis")
  flexsurv_comp("lnorm")
  flexsurv_comp("gengamma")
})