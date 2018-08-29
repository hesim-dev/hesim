context("distributions.R unit tests")
library("flexsurv")

# Methods of moments for beta distribution -------------------------------------
test_that("mom_beta" , {
  beta_mean <- function(a, b) return(a/(a+b))
  beta_var <- function(a, b) return((a * b)/((a + b)^2 * (a + b + 1)))
  beta_params <- mom_beta(.8, .1)
  expect_equal(beta_mean(beta_params$shape1, beta_params$shape2), 
               .8)
  expect_equal(beta_var(beta_params$shape1, beta_params$shape2), 
               .1^2)
  expect_error(mom_beta(.5, sqrt(.5 * (1 - .5))))
})

# Methods of moments for gamma distribution ------------------------------------
test_that("mom_gamma" , {
  # With scale parameter
  gamma_params <- mom_gamma(10000, 1000)
  expect_equal(gamma_params$shape * gamma_params$scale, 
               10000)
  expect_equal(gamma_params$shape * gamma_params$scale^2, 
               1000^2)
  
  # With rate parameter
  gamma_params <- mom_gamma(10000, 1000, scale = FALSE)
  expect_equal(gamma_params$shape / gamma_params$rate,
               10000)
  expect_equal(gamma_params$shape / gamma_params$rate^2,
               1000^2)
})

# Weibull distribution for NMA -------------------------------------------------
shape <-  1.2715
scale <- 6.1914
scalePH <- scale^{-shape}
a0 <- log(shape * scalePH)
a1 <- shape - 1

# pdf
expect_equal(stats::dweibull(.5, shape, scale),
             hesim::dweibullNMA(.5, a0 , a1))
expect_equal(flexsurv::dweibullPH(.5, shape, scalePH),
             hesim::dweibullNMA(.5, a0 , a1))

# cdf
expect_equal(stats::pweibull(4, shape, scale),
             hesim::pweibullNMA(4, a0 , a1))

# quantile
expect_equal(stats::qweibull(.8, shape, scale),
             hesim::qweibullNMA(.8, a0 , a1))
# random
set.seed(3)
r1 <- stats::rweibull(1, shape, scale)
set.seed(3)
r2 <- hesim::rweibullNMA(1, a0, a1)
expect_equal(r1, r2)

# hazard
expect_equal(flexsurv::hweibull(.8, shape, scale),
             hesim::hweibullNMA(.8, a0 , a1))
expect_equal(flexsurv::hweibull(.8, shape, scale, log = TRUE),
             hesim::hweibullNMA(.8, a0 , a1, log = TRUE))

# cumhazard 
expect_equal(flexsurv::Hweibull(.8, shape, scale),
             hesim::HweibullNMA(.8, a0 , a1))
expect_equal(flexsurv::Hweibull(.8, shape, scale, log = TRUE),
             hesim::HweibullNMA(.8, a0 , a1, log = TRUE))

# rmst
expect_equal(flexsurv::rmst_weibull(5, shape, scale),
             hesim::rmst_weibullNMA(5, a0, a1))

# mean
expect_equal(flexsurv::mean_weibull(shape, scale),
             hesim::mean_weibullNMA(a0, a1))

# statistical modeling
## intercept only
s1 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data = bc,
            dist = hesim_survdists$weibullNMA)
s2 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data = bc,
            dist = "weibullPH")
s1.a1 <- s1$res.t["a1", "est"]; s1.a0 <- s1$res.t["a0", "est"]
s2.shape <- s2$res.t["shape", "est"]; s2.scale <- s2$res.t["scale", "est"]
expect_equal(exp(s1.a1), exp(s2.shape) - 1)
expect_equal(s1.a0, log(exp(s2.shape) * exp(s2.scale)))

# covariates on a0
s1 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
            dist = hesim_survdists$weibullNMA)
s2 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
            dist = "weibullPH")
s1.a0 <- s1$res.t["a0", "est"]; s1.a1 <- s1$res.t["a1", "est"]
s2.shape <- s2$res.t["shape", "est"]; s2.scale <- s2$res.t["scale", "est"]
expect_equal(exp(s1.a1), 
             exp(s2.shape) - 1)
expect_equal(log(exp(s2.scale) * exp(s2.shape)), 
             s1.a0)
expect_equal(s2.scale + s2.shape, # since shape is a scalar
             s1.a0)
expect_equal(s1$res.t["groupMedium"],
             s2$res.t["groupMedium"])

# covariates on a0 and a1
s1 <- suppressWarnings(flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
            anc = list(a1 = ~ group),
            dist = hesim_survdists$weibullNMA))
s2 <- suppressWarnings(flexsurvreg(Surv(recyrs, censrec) ~ group, data = bc,
                  anc = list(shape = ~ group),
            dist = "weibullPH"))
s1.a0 <- s1$res.t["a0", "est"]; s1.a1 <- s1$res.t["a1", "est"]
s2.shape <- s2$res.t["shape", "est"]; s2.scale <- s2$res.t["scale", "est"]
expect_equal(exp(s1.a1),
             exp(s2.shape) - 1,
             scale = 1, tol = .001)
expect_equal(s2.scale + s2.shape, # since shape is a scalar
             s1.a0, 
             scale = 1, tol = .001)
expect_equal(exp(sum(s1$res.t[c("a1", "a1(groupMedium)"), "est"])),
             exp(sum(s2$res.t[c("shape", "shape(groupMedium)"), "est"])) - 1,
             scale = 1, tol = .001)

