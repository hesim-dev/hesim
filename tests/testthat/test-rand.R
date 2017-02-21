context("Random number generation")
library("msm")
library("flexsurv")

# Test rcat -------------------------------------------------------------------
test_that("rcat", {
  p <- c(.2, .5, .3, .5, .1, .4)
  n <- 1000
  pmat <- matrix(rep(p, n/2), nrow = n, ncol = 3, byrow = TRUE)
  
  # my function
  set.seed(100)
  samp1 <- rcat(pmat)
  prop.table(table(samp1))
  
  # base R function
  set.seed(100)
  samp2 <- t(apply(pmat, 1, rmultinom, n = 1, size = 1))
  samp2 <- apply(samp2, 1, function(x) which(x == 1))
  prop.table(table(samp2))
  
  # test equal
  expect_equal(samp1, as.matrix(samp2))
})

# Test rpwexp -----------------------------------------------------------------
test_that("rpwexp", {
  rate <- c(.6, 1.2, 1.3)
  n <- 100000
  ratemat <- matrix(rep(rate, n/2), nrow = n, ncol = 3, byrow = TRUE)
  t <- c(0, 10, 15)
  
  # same result as exponential
  set.seed(100)
  samp1 <- rexp(10)
  set.seed(100)
  samp2 <- rpwexp(matrix(rep(1, 10)), t = 0)
  expect_equal(samp1, samp2)
  
  # rpexp
  set.seed(100)
  samp1 <- rpwexp(ratemat, t)
  summary(samp1)
  
  # msm package
  set.seed(100)
  samp2 <- msm::rpexp(n, ratemat[1, ], t)
  summary(samp2) # not identical because of order of sampling from exponential

})

# Test c++ function rsurv -----------------------------------------------------
data(ovarian)
x <- c(1, 55)

test_that("rsurv", {
  n <- 10
  
  ## exponential distribution
  fit <- flexsurv::flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian,
                               dist = "exp")
  fit.lrate <- fit$coef %*% x
  set.seed(50)
  samp1 <- rexp(n, rate = exp(fit.lrate))
  set.seed(50)
  samp2 <- replicate(n, hesim:::rsurv(fit.lrate, anc1 = 0, dist = "exponential"))
  expect_equal(samp1, samp2)
  
  ## weibull distribution
  fit <- flexsurv::flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian,
                               dist = "weibull")
  fit.lscale <- fit$coef[c("scale", "age")] %*% x
  fit.lshape <- fit$coef["shape"]
  set.seed(50)
  samp1 <- rweibull(n, shape = exp(fit.lshape), scale = exp(fit.lscale))
  set.seed(50)
  samp2 <- replicate(n, hesim:::rsurv(fit.lscale, fit.lshape, dist = "weibull"))
  expect_equal(samp1, samp2)
  
  ## gompertz distribution
  fit <- flexsurv::flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian,
                               dist = "gompertz")
  fit.lrate <-  fit$coef[c("rate", "age")] %*% x
  fit.shape <- fit$coef["shape"]
  set.seed(50)
  samp1 <- flexsurv::rgompertz(n, shape = fit.shape, rate = exp(fit.lrate))
  set.seed(50)
  samp2 <- replicate(n, hesim:::rsurv(fit.lrate, fit.shape, dist = "gompertz"))
  expect_equal(samp1, samp2)
  
  ## lognormal distribution
  fit <- flexsurv::flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian,
                               dist = "lnorm")
  fit.meanlog <-  fit$coef[c("meanlog", "age")] %*% x
  fit.lsdlog <- fit$coef["sdlog"]
  set.seed(50)
  samp1 <- rlnorm(n, meanlog = fit.meanlog, sdlog = exp(fit.lsdlog))
  set.seed(50)
  samp2 <- replicate(n, hesim:::rsurv(fit.meanlog, fit.lsdlog, dist = "lnorm"))
  expect_equal(samp1, samp2)
  
  ## gamma distribution
  fit <- flexsurv::flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian,
                               dist = "gamma")
  fit.lrate <- fit$coef[c("rate", "age")] %*% x
  fit.lshape <- fit$coef["shape"]
  set.seed(50)
  samp1 <- rgamma(n, shape = exp(fit.lshape), rate = exp(fit.lrate))
  set.seed(50)
  samp2 <- replicate(n, hesim:::rsurv(fit.lrate, fit.lshape, dist = "gamma"))
  expect_equal(samp1, samp2)
  
  ## log-logistic distribution
  fit <- flexsurv::flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian,
                               dist = "llogis")
  fit.lscale <- fit$coef[c("scale", "age")] %*% x
  fit.lshape <- fit$coef["shape"]
  set.seed(50)
  samp1 <- flexsurv::rllogis(n, shape = exp(fit.lshape), scale = exp(fit.lscale))
  set.seed(50)
  samp2 <- replicate(n, hesim:::rsurv(fit.lscale, fit.lshape, dist = "llogis"))
  expect_equal(samp1, samp2)
  
  ## generalized gamma distribution
  fit <- flexsurv::flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian,
                               dist = "gengamma")
  fit.mu <- fit$coef[c("mu", "age")] %*% x
  fit.lsigma <- fit$coef["sigma"]
  fit.Q <- fit$coef["Q"]
  set.seed(50)
  
  # check that its close to flexsruvreg (not exact because flexsurv function always
  # does lognormal draw)
  samp1 <- flexsurv::rgengamma(n, mu = fit.mu, sigma = exp(fit.lsigma),
                               Q = fit.Q)
  summary(samp1)
  set.seed(50)
  samp2 <- replicate(n, hesim:::rsurv(fit.mu, fit.lsigma, dist = "gengamma", fit.Q))
  summary(samp2)
  
  # lognormal case (Q = 0)
  set.seed(50)
  samp1 <- rlnorm(n, fit.mu, exp(fit.lsigma))
  set.seed(50)
  samp2 <- replicate(n, hesim:::rsurv(fit.mu, fit.lsigma, dist = "gengamma", 0))
  expect_equal(samp1, samp2)
  
  # gamma case (Q = sigma)
  set.seed(50)
  samp1 <- rgamma(n, shape = 1/exp(fit.lsigma)^2, 
                  rate = exp(-fit.mu)/exp(fit.lsigma)^2)
  set.seed(50)
  samp2 <- replicate(n, hesim:::rsurv(fit.mu, fit.lsigma, dist = "gengamma",
                                      exp(fit.lsigma)))
  expect_equal(samp1, samp2)
  
})