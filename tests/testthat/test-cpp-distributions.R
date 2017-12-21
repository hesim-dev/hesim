context("distributions.cpp unit tests")
library("flexsurv")
library("Rcpp")
module <- Rcpp::Module('Distributions', PACKAGE = "hesim")

# Exponential distribution -----------------------------------------------------
test_that("Exponential", {
  Exponential <- module$Exponential
  rate <- 2
  exp <- new(Exponential, rate = rate)
  
  # pdf
  expect_equal(exp$pdf(3), dexp(3, rate = rate))
  
  # cdf 
  expect_equal(exp$cdf(2), pexp(2, rate = rate))
  
  # quantile
  expect_equal(exp$quantile(.025), qexp(.025, rate = rate))
  
  # hazard
  expect_equal(exp$hazard(4), flexsurv::hexp(1, rate = rate))
  
  # cumhazard
  expect_equal(exp$cumhazard(4), flexsurv::Hexp(4, rate = rate))
  
  # random
  set.seed(101)
  r1 <- exp$random()
  set.seed(101)
  r2 <- rexp(1, rate = rate)
  expect_equal(r1, r2)
})

# Weibull distribution ---------------------------------------------------------
test_that("Weibull", {
  Weibull <- module$Weibull
  sh <- 2; sc <- 1.2
  wei <- new(Weibull, shape = sh, scale = sc)

  # pdf
  expect_equal(wei$pdf(3), dweibull(3, shape = sh, scale = sc))

  # cdf
  expect_equal(wei$cdf(2), pweibull(2, shape = sh, scale = sc))

  # quantile
  expect_equal(wei$quantile(.025), qweibull(.025, shape = sh, scale = sc))

  # hazard
  expect_equal(wei$hazard(4), flexsurv::hweibull(4, shape = sh, scale = sc))

  # cumhazard
  expect_equal(wei$cumhazard(4), flexsurv::Hweibull(4, shape = sh, scale = sc))

  # random
  set.seed(101)
  r1 <- wei$random()
  set.seed(101)
  r2 <- rweibull(1, shape = sh, scale = sc)
  expect_equal(r1, r2)
})

# Gamma distribution -----------------------------------------------------------
test_that("Gamma", {
  Gamma <- module$Gamma
  sh <- 2; r <- 1.4
  gamma <- new(Gamma, shape = sh, rate = r)
  
  # pdf
  expect_equal(gamma$pdf(3), dgamma(3, shape = sh, rate = r))
  
  # cdf
  expect_equal(gamma$cdf(2), pgamma(2, shape = sh, rate = r))
  
  # quantile
  expect_equal(gamma$quantile(.025), qgamma(.025, shape = sh, rate = r))
  
  # hazard
  expect_equal(gamma$hazard(4), flexsurv::hgamma(4, shape = sh, rate = r))
  
  # cumhazard
  expect_equal(gamma$cumhazard(4), flexsurv::Hgamma(4, shape = sh, rate = r))
  
  # random
  set.seed(101)
  r1 <- gamma$random()
  set.seed(101)
  r2 <- rgamma(1, shape = sh, rate = r)
  expect_equal(r1, r2)
})

# Lognormal distribution -------------------------------------------------------
test_that("Lognormal", {
  Lognormal <- module$Lognormal
  m <- 8; s <- 2.5
  lognormal <- new(Lognormal, meanlog = m, sdlog = s)
  
  # pdf
  expect_equal(lognormal$pdf(3), dlnorm(3, meanlog = m, sdlog = s))
  
  # cdf
  expect_equal(lognormal$cdf(2), plnorm(2, meanlog = m, sdlog = s))
  
  # quantile
  expect_equal(lognormal$quantile(.025), qlnorm(.025, meanlog = m, sdlog = s))
  
  # hazard
  expect_equal(lognormal$hazard(4), flexsurv::hlnorm(4, meanlog = m, sdlog = s))
  
  # cumhazard
  expect_equal(lognormal$cumhazard(4), flexsurv::Hlnorm(4, meanlog = m, sdlog = s))
  
  # random
  set.seed(101)
  r1 <- lognormal$random()
  set.seed(101)
  r2 <- rlnorm(1, meanlog = m, sdlog = s)
  expect_equal(r1, r2)
})

# Gompertz distribution --------------------------------------------------------
test_that("Gompertz", {
  Gompertz <- module$Gompertz
  
  # shape > 0
  sh <- .05; r <- .5
  gompertz <- new(Gompertz, shape = sh, rate = r)
  
  ## pdf
  expect_equal(gompertz$pdf(2), flexsurv::dgompertz(2, shape = sh, rate = r))
  
  ## cdf
  expect_equal(gompertz$cdf(2), flexsurv::pgompertz(2, shape = sh, rate = r))
  
  ## quantile
  expect_equal(gompertz$quantile(.35), flexsurv::qgompertz(.35, shape = sh, rate = r))

  ## hazard
  expect_equal(gompertz$hazard(4), flexsurv::hgompertz(4, shape = sh, rate = r))
  
  ## cumhazard
  expect_equal(gompertz$cumhazard(4), flexsurv::Hgompertz(4, shape = sh, rate = r))
  
  ## random
  set.seed(101)
  r1 <- gompertz$random()
  set.seed(101)
  r2 <- flexsurv::rgompertz(1, shape = sh, rate = r)
  expect_equal(r1, r2)
  
  # shape = 0
  sh <- 0; r <- .5
  gompertz <- new(Gompertz, shape = sh, rate = r)
  
  ## pdf
  expect_equal(gompertz$pdf(2), flexsurv::dgompertz(2, shape = sh, rate = r))
  
  ## cdf
  expect_equal(gompertz$cdf(2), flexsurv::pgompertz(2, shape = sh, rate = r))
  
  # quantile
  expect_equal(gompertz$quantile(.7), flexsurv::qgompertz(.7, shape = sh, rate = r))
  
  # hazard
  expect_equal(gompertz$hazard(.7), flexsurv::hgompertz(.7, shape = sh, rate = r))
  
  # cumhazard
  expect_equal(gompertz$cumhazard(.7), flexsurv::Hgompertz(.7, shape = sh, rate = r))
})

# Log-logistic distribution ----------------------------------------------------
test_that("LogLogistic", {
  LogLogistic <- module$LogLogistic
  sh <- 1; sc <- .5
  llogis <- new(LogLogistic, shape = sh, scale = sc)
  
  # pdf
  expect_equal(llogis$pdf(2), flexsurv::dllogis(2, shape = sh, scale = sc))
  
  # cdf
  expect_equal(llogis$cdf(2), flexsurv::pllogis(2, shape = sh, scale = sc))

  # quantile
  expect_equal(llogis$quantile(.34), flexsurv::qllogis(.34, shape = sh, scale = sc))
  
  # hazard
  expect_equal(llogis$hazard(6), flexsurv::hllogis(6, shape = sh, scale = sc))
  
  # cumhazard
  expect_equal(llogis$cumhazard(6), flexsurv::Hllogis(6, shape = sh, scale = sc))
  
  # random
  set.seed(101)
  r1 <- llogis$random()
  set.seed(101)
  r2 <- flexsurv::rllogis(1, shape = sh, scale = sc)
  expect_equal(r1, r2)
  
})

# Generalized gamma distribution -----------------------------------------------
test_that("GeneralizedGamma", {
  GeneralizedGamma <- module$GeneralizedGamma
  
  # Q < 0
  m <- 2; s <- 1.5; q <- -2
  gengamma <- new(GeneralizedGamma, mu = m, sigma = s, Q = q)
  
  ## pdf
  expect_equal(gengamma$pdf(4), flexsurv::dgengamma(4, mu = m, sigma = s, Q = q))
  
  ## cdf
  expect_equal(gengamma$cdf(4), flexsurv::pgengamma(4, mu = m, sigma = s, Q = q))
  
  ## quantile
  expect_equal(gengamma$quantile(.5), flexsurv::qgengamma(.5, mu = m, sigma = s, Q = q))
  
  ## hazard
  expect_equal(gengamma$hazard(1), flexsurv::hgengamma(1, mu = m, sigma = s, Q = q))
  
  ## cumhazard
  expect_equal(gengamma$cumhazard(1), flexsurv::Hgengamma(1, mu = m, sigma = s, Q = q))
  
  # Q = 0
  m <- 5; s <- 2; q <- 0
  gengamma <- new(GeneralizedGamma, mu = m, sigma = s, Q = q)
  
  ## pdf
  expect_equal(gengamma$pdf(4), flexsurv::dgengamma(4, mu = m, sigma = s, Q = q))
  
  ## cdf
  expect_equal(gengamma$cdf(4), flexsurv::pgengamma(4, mu = m, sigma = s, Q = q))
  
  ## quantile
  expect_equal(gengamma$quantile(.8), flexsurv::qgengamma(.8, mu = m, sigma = s, Q = q))
  
  ## hazard
  expect_equal(gengamma$hazard(1), flexsurv::hgengamma(1, mu = m, sigma = s, Q = q))
  
  ## cumhazard
  expect_equal(gengamma$cumhazard(1), flexsurv::Hgengamma(1, mu = m, sigma = s, Q = q))
  
  ## random
  set.seed(101)
  r1 <- gengamma$random()
  set.seed(101)
  r2 <- flexsurv::rgengamma(1, mu = m, sigma = s, Q = q)
  expect_equal(r1, r2)
  
  # Q > 0
  m <- 2; s <- 1.1; q <- 2
  gengamma <- new(GeneralizedGamma, mu = m, sigma = s, Q = q)
  
  ## pdf
  expect_equal(gengamma$pdf(2), flexsurv::dgengamma(2, mu = m, sigma = s, Q = q))
  
  ## cdf
  expect_equal(gengamma$cdf(4), flexsurv::pgengamma(4, mu = m, sigma = s, Q = q))
  
  ## quantile
  expect_equal(gengamma$quantile(.8), flexsurv::qgengamma(.8, mu = m, sigma = s, Q = q))
  
  ## hazard
  expect_equal(gengamma$hazard(1), flexsurv::hgengamma(1, mu = m, sigma = s, Q = q))
  
  ## cumhazard
  expect_equal(gengamma$cumhazard(1), flexsurv::Hgengamma(1, mu = m, sigma = s, Q = q))

})


# Truncated normal distribution ------------------------------------------------
test_that("rtruncnorm", {
  n <- 1000
  mu <- 50; sigma <- 10; lower <- 25; upper <- 60
  
  #rtruncnorm from hesim
  set.seed(10)
  samp1 <- replicate(n, hesim:::rtruncnorm(mu, sigma, lower, upper))

  # rtruncnorm from truncnorm package
  set.seed(10)
  samp3 <- truncnorm::rtruncnorm(n, lower, upper, mu, sigma)
  expect_equal(samp1, samp3)
})