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
  expect_equal(length(fit.sample), length(fit$dlist$pars))
  expect_equal(c(fit.sample$shape), fit.normboot[, "shape"])
  expect_equal(nrow(fit.sample$scale), n)
  expect_false(any(lapply(fit.sample[1:2], ncol) != 1))
  expect_true(inherits(fit.sample, "surv_pars"))
  
  ## sample of size 1
  expect_error(rsample(fit, n = 1), NA)
  
  # covariates on 1 paramters
  ## lognormal
  fit <- flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                     data = ovarian, dist = "lognormal")
  fit.sample <- rsample(fit, n = 3)
  expect_equal(ncol(fit.sample$meanlog), 2)
  expect_equal(ncol(fit.sample$sdlog), 1)
  
  # covariates on 2 parameters
  fit <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                     anc = list(sigma = ~ group), dist = "gengamma") 
  fit.sample <- rsample(fit, n = 2)
  expect_equal(ncol(fit.sample$mu), 3)
  expect_equal(ncol(fit.sample$sigma), 3)
  expect_equal(ncol(fit.sample$Q), 1)
  
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
  expect_equal(names(fit2.samp), names(dm))
  
  # list
  fit.list <- list(fit1, fit2)
  dm <- design_matrix(fit.list)
  expect_equal(length(dm), 2)
  expect_true(names(dm[[1]]) == "rate")
  expect_equal(dm[[2]]$rate[1, "age"], fit2$data$m$age[1])
})

# Survival R6 class ------------------------------------------------------------
test_that("DisModSurv class no regressors", {
  fit <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                     data = ovarian, dist = "exponential")
  surv <- DisModSurv$new(dist_name = fit$dlist$name,
                       coefs = rsample(fit, n = 5),
                       X = design_matrix(fit))

  # Quantiles
  surv.quantiles <- surv$quantiles(q = c(.5, .8))
  expect_equal(stats::qexp(.5, exp(surv$coefs$rate[1])), 
               surv.quantiles[sim == 1 & id == 1 & q == .5, value])
  
  # Survival curves
  surv.sc <- surv$surv(t = c(0, 1, 2))
  expect_equal(1 - stats::pexp(1, rate = exp(surv$coefs$rate[2])),
                surv.sc[sim == 2 & id == 1 & t == 1, value])
  
  # Cummulative hazards
  surv.cumhaz <- surv$cumhazard(t = c(0, 1, 2))
  expect_equal(flexsurv::Hexp(2, rate = exp(surv$coefs$rate[3])),
               surv.cumhaz[sim == 3 & id == 1 & t == 2, value])
  
  # Hazards
  surv.haz <- surv$hazard(t = c(0, 1, 2, 3))
  expect_equal(flexsurv::hexp(3, rate = exp(surv$coefs$rate[3])),
               surv.haz[sim == 3 & id == 1 & t == 3, value])
  
  # rmst
  expect_true(surv$rmst(t = c(15, 20), r = 0)[1, "value"] > 
                surv$rmst(t = c(15, 20), r = 0.03)[1, "value"])
  
  # set inputs
  surv$set_inputs(X = 1, coefs = 1)
  expect_equal(surv$X, 1)
  expect_equal(surv$coefs, 1)
  expect_error(surv$set_inputs(param = 3))
  expect_error(surv$set_inputs(coefs = 3, param = 3))
  expect_error(surv$set_inputs(rmst = 3))
  expect_error(surv$set_inputs(rmst = function(x) x))
})

test_that("DisModSurv class with regressors", {
  fit <- flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                     data = ovarian, dist = "weibull")
  surv <- DisModSurv$new(dist_name = fit$dlist$name,
                       coefs = rsample(fit, n = 5),
                       X = design_matrix(fit))
  
  # Quantiles
  surv.quantiles <- surv$quantiles(q = c(.5))
  scale1 <- surv$X$scale[1,, drop = FALSE] %*% t(surv$coefs$scale[1,, drop = FALSE])
  expect_equal(stats::qweibull(.5, shape = exp(surv$coefs$shape[1]), 
                               scale = exp(scale1)), 
               surv.quantiles[sim == 1 & id == 1 & q == .5, value])
  
  # Survival curves
  surv.sc <- surv$surv(t = c(100, 365))
  expect_equal(1 - stats::pweibull(365, shape = exp(surv$coefs$shape[1]), 
                               scale = exp(scale1)), 
               surv.sc[sim == 1 & id == 1 & t == 365, value])
})

test_that("Compare survival class with flexsurv", {
  flexsurv_comp <- function(dist_name){
      fit <- flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                        data = ovarian, dist = dist_name)
      surv <- DisModSurv$new(dist_name = fit$dlist$name,
                           coefs = rsample(fit, point_estimate = TRUE),
                            X = design_matrix(fit, indices = 1))
    
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

# Survival Model R6 class ------------------------------------------------------
test_that("DecModSurv class no regressors", {
  fit <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
                     data = ovarian, dist = "weibull")
  dis.surv <- DisModSurv$new(dist_name = fit$dlist$name,
                      coefs = rsample(fit, n = 5),
                      X = design_matrix(fit),
                      time_length = 1/365)
  util.vals <- UtilityValues$new(mean = rep(1, 5))
  dec.surv <- DecModSurv$new(dis_mod_surv = dis.surv,
                                utility_values = util.vals)

  # Simulate Effects
  ## utility = 1
  dec.surv$sim_effects(max_t = 700, r = c(0, .03), type = c("lys", "qalys"))
  dis.surv.rmst <- dis.surv$rmst(t = 700, r = 0)
  expect_equal(as.numeric(dec.surv$effects_[sim == 3 & id == 2 & type == "lys", "value"]),
               as.numeric(dis.surv.rmst[sim == 3 & id == 2, "value"]))
  
  dec.surv$sim_effects(max_t = 700, r = .03, type = "lys")
  dis.surv.rmst <- dis.surv$rmst(t = 700, r = .03)
  expect_equal(as.numeric(dec.surv$effects_[sim == 2 & id == 3 & type == "lys", "value"]),
               as.numeric(dis.surv.rmst[sim == 2 & id == 3, "value"]))
  
  # utility = .8
  dec.surv$utility_values <- UtilityValues$new(mean = rep(.8, 5))
  dec.surv$sim_effects(max_t = 700, r = .03, type = "qalys")
  test_integrate <- function(w, r = .03, upper){
    f <- function(t) {
      w * exp(-r * t/365) * (1 - pweibull(t, shape = exp(dis.surv$coefs$shape[1]),
                                           scale = exp(dis.surv$coefs$scale[1])))
    }
    stats::integrate(f, 0, upper)$value
  }
  expect_equal(as.numeric(dec.surv$effects_[sim == 1 & id == 1 & type == "qalys", "value"]),
               test_integrate(.8, r = .03, upper = 700), 
               tolerance = .001, scale = 1)
  
  # Simulate costs
  dec.surv$cost_values <- CostValues$new(list(matrix(rep(10, 5)),
                                              matrix(rep(12, 5))),
                                         names = c("cost1", "cost2"))
  dec.surv$sim_costs(t = c(14, 22), r = .02)
  expect_equal(dec.surv$costs_[sim == 1 & id == 1 & component == "cost1"]$value,
                 test_integrate(10, r = .02, upper = 14),
                   tolerance = .001, scale = 1)
  expect_equal(dec.surv$costs_[sim == 1 & id == 1 & component == "cost2"]$value,
               test_integrate(12, r = .02, upper = 22),
               tolerance = .001, scale = 1)
  
  # set inputs
  dec.surv$set_inputs(cost_values = 20, utility_values = 2)
  expect_equal(dec.surv$cost_values, 20)
  expect_equal(dec.surv$utility_values, 2)
  expect_error(dec.surv$set_inputs(param = 3))
  expect_error(dec.surv$set_inputs(utility_values = 3, param = 3))
})

test_that("DecModSurv class with regressors", {
  fit <- flexsurvreg(formula = Surv(futime, fustat) ~ age, 
                     data = ovarian, dist = "weibull")
  dis.surv <- DisModSurv$new(dist_name = fit$dlist$name,
                      coefs = rsample(fit, n = 5),
                      X = design_matrix(fit),
                      time_length = 1/365)
  util.vals <- UtilityValues$new(mean = rep(1, 5))
  dec.surv <- DecModSurv$new(dis_mod_surv = dis.surv,
                                utility_values = util.vals)
  
  # Simulate Effects
  dec.surv$sim_effects(max_t = 700, r = c(0, .03), type = c("lys", "qalys"))
  dis.surv.rmst <- dis.surv$rmst(t = 700, r = 0)
  expect_equal(as.numeric(dec.surv$effects_[sim == 3 & id == 2 & type == "lys", "value"]),
               as.numeric(dis.surv.rmst[sim == 3 & id == 2, "value"]))
  
  dec.surv$sim_effects(max_t = 700, r = .03, type = "lys")
  dis.surv.rmst <- dis.surv$rmst(t = 700, r = .03)
  expect_equal(as.numeric(dec.surv$effects_[sim == 2 & id == 3 & type == "lys", "value"]),
               as.numeric(dis.surv.rmst[sim == 2 & id == 3, "value"]))
  
  # Simulate Costs
  # surv.mod$sim_costs(t = 700, r = 0)$costs_
})
