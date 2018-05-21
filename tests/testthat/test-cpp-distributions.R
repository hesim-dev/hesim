context("distributions.cpp unit tests")
library("flexsurv")
library("numDeriv")
library("Rcpp")
module <- Rcpp::Module('Distributions', PACKAGE = "hesim")

# Free functions ---------------------------------------------------------------
shape <- 2
scale <- 3
parsC <- hesim:::C_weibull_to_weibullNMA(shape, scale)
parsR <- hesim:::weibull_to_weibullNMA(shape, scale)
expect_equal(parsC[1], parsR$a0)
expect_equal(parsC[2], parsR$a1)

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

# Survival splines -------------------------------------------------------------
basis_cube <- function(x){
  return (max(0, x^3))
}

R_linear_predict <- function(t, gamma, knots, timescale){
  res <- rep(NA, length(t))
  for (k in 1:length(t)){
    t.scaled <- switch(timescale,
                     log = log(t[k]),
                     identity = t[k])
    knot.min <- knots[1];  knot.max <- knots[length(knots)]
    basis <- rep(NA, length(knots))
    basis[1] <- 1; basis[2] <- t.scaled
    for (j in 2:(length(knots) - 1)){
      lambda.j <- (knot.max - knots[j])/(knot.max - knot.min)
      basis[j + 1] <- basis_cube(t.scaled - knots[j]) - 
                      lambda.j * basis_cube(t.scaled - knot.min) -
                      (1 - lambda.j) * basis_cube(t.scaled - knot.max)
    }
    res[k] <- basis %*% gamma
  }
  return (res)
}

test_that("SurvivalSplines", {
  SurvSplines <- module$SurvSplines

  # log hazard
  ## Test function
  test_SurvSplines1 <- function(surv_splines, timescale, 
                                gamma = c(-1.2, 1.3, .07),
                                knots = c(.19, 1.7, 6.7)){ 
    surv.splines <- new(surv_splines, gamma = gamma, knots = knots,
                          scale = "log_hazard", timescale = timescale)
    t <- 2
    linear.predict <- surv.splines$linear_predict(t)
    expect_equal(linear.predict,
                 R_linear_predict(t, gamma = gamma, knots = knots, 
                                  timescale = timescale))
    linear.predict.dx <- surv.splines$linear_predict_dx(t)
      
    ### hazard
    expect_equal(surv.splines$hazard(t), 
                   exp(linear.predict))
    
    ### cummulative hazard
    R_hazard <- function(t) {
      exp(R_linear_predict(t, gamma = gamma, knots = knots,
                                          timescale = timescale))
    }
    R_cumhazard <- function(t){
      stats::integrate(R_hazard, 0, t)$value
    }
    expect_equal(surv.splines$cumhazard(t), R_cumhazard(t),
                 tolerance = .001, scale = 1)
    
    ### cdf
    expect_equal(surv.splines$cdf(.8),
                 1 - exp(-surv.splines$cumhazard(.8)))
    
    ### pdf
    R_cdf <- function(t){
      return(1 - exp(-R_cumhazard(t)))
    }
    expect_equal(surv.splines$pdf(0), 0)
    expect_equal(surv.splines$pdf(.5), 
                 (1 - surv.splines$cdf(.5)) * surv.splines$hazard(.5))
    expect_equal(surv.splines$pdf(.5),
                 numDeriv::grad(R_cdf, .5))
    
    ### quantile
    expect_equal(surv.splines$cdf(surv.splines$quantile(.55)), .55,
                 tolerance = .001, scale = 1)
    
  }
  test_SurvSplines1(SurvSplines, timescale = "identity")
  test_SurvSplines1(SurvSplines, timescale = "log")
  
  # log cummulative hazard, log odds, and normal
  ## test function
  test_SurvSplines2 <- function(surv_splines, flexsurv_scale, timescale,
                                gamma = NULL, knots = NULL){
    if (is.null(gamma) | is.null(knots)){
      spl.fit <- flexsurvspline(Surv(recyrs, censrec) ~ 1, 
                          data = bc, k = 1, timescale = timescale,
                          scale = flexsurv_scale)
      gamma <- spl.fit$res.t[, "est"]
      knots <- spl.fit$knots
    }
    hesim.scale <- switch(flexsurv_scale,
                          hazard = "log_cumhazard",
                          odds = "log_cumodds",
                          normal = "inv_normal"
                          )
    
    surv.splines <- new(surv_splines, gamma = gamma, knots = knots,
                    scale = hesim.scale, timescale = timescale)
    
    ### pdf
    expect_equal(dsurvspline(5, gamma = gamma, knots = knots, 
                             scale = flexsurv_scale, timescale = timescale),
                surv.splines$pdf(5))
    expect_equal(surv.splines$pdf(0), 0)
    expect_equal(surv.splines$pdf(-5), 0)
    
    ### cdf
    expect_equal(psurvspline(.5, gamma = gamma, knots = knots, 
                             scale = flexsurv_scale, timescale = timescale),
               surv.splines$cdf(.5))
    expect_equal(surv.splines$cdf(0), 0)
    expect_equal(surv.splines$cdf(-5), 0)
    
    ### quantile
    expect_equal(qsurvspline(.8, gamma = gamma, knots = knots, 
                             scale = flexsurv_scale, timescale = timescale),
                  surv.splines$quantile(.8), tolerance = .001, scale = 1)
    expect_equal(surv.splines$quantile(-2), NaN)
    expect_equal(surv.splines$quantile(0), -Inf)
    expect_equal(surv.splines$quantile(1), Inf)
    
    ### hazard
    expect_equal(hsurvspline(3, gamma = gamma, knots = knots, 
                             scale = flexsurv_scale, timescale = timescale),
                  surv.splines$hazard(3))
    expect_equal(surv.splines$hazard(0), 0)
    expect_equal(surv.splines$hazard(-5), 0)
    
    ### cumhazard
    expect_equal(Hsurvspline(3, gamma = gamma, knots = knots, 
                             scale = flexsurv_scale, timescale = timescale),
                 surv.splines$cumhazard(3))
    expect_equal(surv.splines$cumhazard(0), 0)
    expect_equal(surv.splines$cumhazard(-5), 0)
    
    ### random
    set.seed(12)
    r1 <- rsurvspline(1, gamma = gamma, knots = knots, 
                      scale = flexsurv_scale, timescale = timescale)
    set.seed(12)
    r2 <- surv.splines$random()
    expect_equal(r1, r2, tolerance = .001, scale = 1)
  }
  test_SurvSplines2(surv_splines = SurvSplines, flexsurv_scale = "hazard",
                    timescale = "log")
  test_SurvSplines2(surv_splines = SurvSplines, flexsurv_scale = "hazard",
                    timescale = "identity")
  test_SurvSplines2(surv_splines = SurvSplines, flexsurv_scale = "odds",
                    timescale = "log")
  test_SurvSplines2(surv_splines = SurvSplines, flexsurv_scale = "odds",
                    timescale = "identity")
  test_SurvSplines2(surv_splines = SurvSplines, flexsurv_scale = "normal",
                    timescale = "log")
  test_SurvSplines2(surv_splines = SurvSplines, flexsurv_scale = "normal",
                    timescale = "identity", gamma = c(-1.2, 1.3, .07),
                    knots = c(.19, 1.7, 6.7))
})

# Survival Fractional Polynomials ----------------------------------------------
# functions from flexsurv tests
bfp <- function (x, powers = c(1, 2)) {
  nobs <- length(x)
  npoly <- length(powers)
  X <- matrix(0, nrow = nobs, ncol = npoly)
  x1 <- ifelse(powers[1] != rep(0, nobs), x^powers[1], log(x))
  X[, 1] <- x1
  if (npoly >= 2) {
      for (i in 2:npoly) {
          if (powers[i] == powers[(i - 1)]) 
              x2 <- log(x) * x1
          else x2 <- ifelse(powers[i] != rep(0, nobs), x^powers[i], 
              log(x))
          X[, i] <- x2
          x1 <- x2
      }
  }
  X
}

hfp.lh <- function(x, gamma, powers){
  if(!is.matrix(gamma)) gamma <- matrix(gamma, nrow=1)
  lg <- nrow(gamma)
  nret <- max(length(x), lg)
  gamma <- apply(gamma, 2, function(x)rep(x,length=nret))
  x <- rep(x, length=nret)
  basis <- cbind(1, bfp(x, powers))
  loghaz <- rowSums(basis * gamma)
  exp(loghaz)
}

hfp.lh3 <- unroll.function(hfp.lh, gamma = 0:2)

custom.hfp.lh3 <- list(
  name = "fp.lh3",
  pars = c(paste0("gamma", 0:2)),
  location = c("gamma0"),
  transforms = rep(c(identity), 3), inv.transforms = rep(c(identity), 3)
)

# fit fractional polynomial model
powers <- c(1, 0)
# fp.fit <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data = bc, 
#                       aux = list(powers = powers), inits = c(-2, 0, 0),
#                       dist = custom.hfp.lh3)
gamma <- c(-1.2, -.567, 1.15)

test_that("FracPoly", {
  FracPoly <- module$FracPoly
  
  # specifications 1
  fp <- new(FracPoly, gamma = gamma, powers = powers)
  R_hazard <- function(t){
    return(exp(c(gamma %*% t(cbind(1, bfp(t, powers))))))
  }
  R_cumhazard <- function(t){
    stats::integrate(R_hazard, 0, t)$value
  }
  
  ## linear prediction
  expect_equal(fp$linear_predict(3), 
               c(gamma %*% t(cbind(1, bfp(3, powers)))))
  
  ## hazard
  expect_equal(fp$hazard(4), 
               exp(fp$linear_predict(4)))
  expect_equal(fp$hazard(0), 0)
  
  ## cumhazard
  expect_equal(fp$cumhazard(2), 
               R_cumhazard(2),
               tolerance = .001, scale = 1)
  expect_equal(fp$cumhazard(0), 0)
  
  ## cdf
  expect_equal(fp$cdf(2.5), 
               1 - exp(-fp$cumhazard(2.5)))
  
  ## pdf
  R_cdf <- function(t){
    1 - exp(-R_cumhazard(t))
  }
  expect_equal(fp$pdf(.8), 
               numDeriv::grad(R_cdf, .8))
  
  ## quantile
  expect_equal(.45, 
               fp$cdf(fp$quantile(.45)),
               tol = .001, scale = 1)
  
  ## random
  set.seed(101)
  r1 <- fp$random()
  set.seed(101)
  r2 <- fp$quantile(runif(1, 0, 1))
  expect_equal(r1, r2)
  
  # other specifications
  test_linear_predict <- function(gamma, powers){
    fp <- new(FracPoly, gamma = gamma, powers = powers)
    expect_equal(fp$linear_predict(3), 
                        c(gamma %*% t(cbind(1, bfp(3, powers)))))
  }
  test_linear_predict(c(.2, .2), 0)
  test_linear_predict(c(.2, .2, .2), c(0, 0))
  test_linear_predict(c(.2, .2, .2), c(1, 1))
  test_linear_predict(c(.2, .2, .5, .5, .2, .2), c(1, 1, 2, 2, 1))
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
