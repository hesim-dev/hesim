context("distributions.cpp unit tests")
library("flexsurv")
library("numDeriv")
library("Rcpp")
module <- Rcpp::Module('distributions', PACKAGE = "hesim")

# Exponential distribution -----------------------------------------------------
test_that("exponential", {
  Exponential <- module$exponential
  rate <- 2
  exp <- new(Exponential, rate = rate)
  
  # pdf
  expect_equal(exp$pdf(3), 
               dexp(3, rate = rate))
  
  # cdf
  expect_equal(exp$cdf(3), 
               pexp(3, rate = rate))  

  # quantile
  expect_equal(exp$quantile(.025), 
               qexp(.025, rate = rate))    
  
  # hazard
  expect_equal(exp$hazard(1), 
               flexsurv::hexp(1, rate = rate))  
  
  # cumhazard
  expect_equal(exp$cumhazard(4), 
               flexsurv::Hexp(4, rate = rate))  
  
  # random
  set.seed(101)
  r1 <- exp$random()
  set.seed(101)
  r2 <- rexp(1, rate = rate)
  expect_equal(r1, r2)
  
  # Truncated random
  r <- replicate(10, exp$trandom(1, 5))
  expect_true(all(r >= 1))
  expect_true(all(r <= 5))
})

# Piecewise exponential distribution -------------------------------------------
test_that("pwexp", {
  PwExp <- module$piecewise_exponential
  rate <- c(.8, 1.2, 4)
  time <- c(0, 5, 9)
  pwexp <- new(PwExp, rate = rate, time = time)
  
  # pdf
  expect_equal(pwexp$pdf(5.5), 1.2*exp(-(5.0*0.8 + 0.5*1.2)))
  
  # cdf
  expect_equal(pwexp$cdf(5.5), 1.0-exp(-(5.0*0.8 + 0.5*1.2)))
  
  # quantile
  expect_equal(pwexp$quantile(1.0-exp(-(5.0*0.8 + 0.5*1.2))), 5.5)
  
  # hazard
  expect_equal(pwexp$hazard(5.5), 1.2)
  expect_equal(pwexp$hazard(0.0), 0.8)
  
  # cumhazard
  expect_equal(pwexp$cumhazard(5.5), 5.0*0.8 + 0.5*1.2)
  expect_equal(pwexp$cumhazard(5.0), 5.0*0.8)
  
  # random
  set.seed(101)
  r1 <- pwexp$random()
  set.seed(101)
  r2 <- rpwexp(1, rate = rate, time = time)
  expect_equal(r1, r2)
  
  # Truncated random
  expect_error(
    pwexp$trandom(1, 5),
    paste0("hesim does not currently support sampling from a piecewise ", 
           "exponential distribution truncated from above.")
  )
  
  r <- replicate(10, pwexp$trandom(7, Inf))
  expect_true(all(r >= 7))
  
  r <- replicate(10, pwexp$trandom(max(time) + 1, Inf))
  expect_true(all(r >= max(time) + 1))
})

# Weibull distribution (AFT) ---------------------------------------------------
test_that("weibull", {
  Weibull <- module$weibull
  sh <- 2; sc <- 1.2
  wei <- new(Weibull, shape = sh, scale = sc)
  
  # pdf
  expect_equal(wei$pdf(3), 
               dweibull(3, shape = sh, scale = sc)) 
  
  # cdf
  expect_equal(wei$cdf(2), 
               pweibull(2, shape = sh, scale = sc))   
  
  # quantile
  expect_equal(wei$quantile(.025), 
               qweibull(.025, shape = sh, scale = sc))    
  
  # hazard
  expect_equal(wei$hazard(4), 
               flexsurv::hweibull(4, shape = sh, scale = sc))    
  
  # cumhazard
  expect_equal(wei$cumhazard(4), 
               flexsurv::Hweibull(4, shape = sh, scale = sc))  
  
  # random
  set.seed(101)
  r1 <- wei$random()
  set.seed(101)
  r2 <- rweibull(1, shape = sh, scale = sc)
  expect_equal(r1, r2)
  
  # Truncated random
  r <- replicate(10, wei$trandom(0, 1))
  expect_true(all(r >= 0))
  expect_true(all(r <= 1))  
})

# Weibull distribution (PH) ----------------------------------------------------
test_that("weibull_ph", {
  WeibullPH <- module$weibull_ph
  a <- exp(0.36); m <- exp(-5)
  wei <- new(WeibullPH, shape = a, scale = m)
  
  # pdf
  expect_equal(wei$pdf(4), 
               flexsurv::dweibullPH(4, shape = a, scale = m))   
  
  # cdf
  expect_equal(wei$cdf(4), 
               flexsurv::pweibullPH(4, shape = a, scale = m))   
  
  # quantile
  expect_equal(wei$quantile(.7), 
               flexsurv::qweibullPH(.7, shape = a, scale = m))     
  
  # hazard
  expect_equal(wei$hazard(2), 
               flexsurv::hweibullPH(2, shape = a, scale = m))     
  
  # cumhazard
  expect_equal(wei$cumhazard(2), 
               flexsurv::HweibullPH(2, shape = a, scale = m))       
  
  # random
  set.seed(101)
  r1 <- wei$random()
  set.seed(101)
  r2 <- rweibullPH(1, shape = a, scale = m)
  expect_equal(r1, r2)
  
  # Truncated random
  r <- replicate(10, wei$trandom(2, 11))
  expect_true(all(r >= 2))
  expect_true(all(r <= 11))  
})

# Weibull distribution for NMA -------------------------------------------------
test_that("weibull_nma", {
  WeibullNma <- module$weibull_nma
  a0 <- -2.07; a1 <- .2715
  wei <- new(WeibullNma, a0 = a0, a1 = a1)
  
  # pdf
  expect_equal(wei$pdf(4), 
               dweibullNMA(4, a0 = a0, a1 = a1))   
  
  # cdf
  expect_equal(wei$cdf(4), 
               pweibullNMA(4, a0 = a0, a1 = a1))   
  
  # quantile
  expect_equal(wei$quantile(.7), 
               qweibullNMA(.7, a0 = a0, a1 = a1))     
  
  # hazard
  expect_equal(wei$hazard(2), 
               hweibullNMA(2, a0 = a0, a1 = a1))     
  
  # cumhazard
  expect_equal(wei$cumhazard(2), 
               HweibullNMA(2, a0 = a0, a1 = a1))       
  
  # random
  set.seed(101)
  r1 <- wei$random()
  set.seed(101)
  r2 <- rweibullNMA(1, a0 = a0, a1 = a1)
  expect_equal(r1, r2)
  
  # Truncated random
  r <- replicate(10, wei$trandom(0, 11))
  expect_true(all(r >= 0))
  expect_true(all(r <= 11))  
})

# Gamma distribution -----------------------------------------------------------
test_that("gamma", {
  Gamma <- module$gamma
  sh <- 2; r <- 1.4
  gamma <- new(Gamma, shape = sh, rate = r)
  
  # pdf
  expect_equal(gamma$pdf(3), 
               dgamma(3, shape = sh, rate = r))   
  
  # cdf
  expect_equal(gamma$cdf(2), 
               pgamma(2, shape = sh, rate = r))     
  
  # quantile
  expect_equal(gamma$quantile(.025), 
               qgamma(.025, shape = sh, rate = r))       
  
  # hazard
  expect_equal(gamma$hazard(4), 
               flexsurv::hgamma(4, shape = sh, rate = r))     
  
  # cumhazard
  expect_equal(gamma$cumhazard(4), 
               flexsurv::Hgamma(4, shape = sh, rate = r))      
  
  # random
  set.seed(101)
  r1 <- gamma$random()
  set.seed(101)
  r2 <- rgamma(1, shape = sh, rate = r)
  expect_equal(r1, r2)
  
  # Truncated random
  r <- replicate(10, gamma$trandom(0, 11))
  expect_true(all(r >= 0))
  expect_true(all(r <= 11))   
})

# Lognormal distribution -------------------------------------------------------
test_that("lognormal", {
  Lognormal <- module$lognormal
  m <- 8; s <- 2.5
  lnorm <- new(Lognormal, meanlog = m, sdlog = s)
  
  # pdf
  expect_equal(lnorm$pdf(3), 
               dlnorm(3, meanlog = m, sdlog = s))   
  
  # cdf
  expect_equal(lnorm$cdf(3), 
               plnorm(3, meanlog = m, sdlog = s))     
  
  # quantile
  expect_equal(lnorm$quantile(.33), 
               qlnorm(.33, meanlog = m, sdlog = s))     
  
  # hazard
  expect_equal(lnorm$hazard(4), 
               flexsurv::hlnorm(4, meanlog = m, sdlog = s))    
  
  # cumhazard
  expect_equal(lnorm$cumhazard(4), 
               flexsurv::Hlnorm(4, meanlog = m, sdlog = s))     
  
  # random
  set.seed(101)
  r1 <- lnorm$random()
  set.seed(101)
  r2 <- rlnorm(1, meanlog = m, sdlog = s)
  expect_equal(r1, r2)
  
  # Truncated random
  r <- replicate(10, lnorm$trandom(10, 21))
  expect_true(all(r >= 10))
  expect_true(all(r <= 21))     
})

# Gompertz distribution --------------------------------------------------------
test_that("gompertz", {
  Gompertz <- module$gompertz  
  
  test_gompertz <- function(sh, r){
    
    gomp <- new(Gompertz, shape = sh, rate = r)
      
    ## pdf
    expect_equal(gomp$pdf(2), 
                 flexsurv::dgompertz(2, shape = sh, rate = r))   
    
    ## cdf
    expect_equal(gomp$cdf(1.5), 
                 flexsurv::pgompertz(1.5, shape = sh, rate = r))   
  
    ## quantile
    expect_equal(gomp$quantile(.21), 
                 flexsurv::qgompertz(.21, shape = sh, rate = r))     
    
    ## hazard
    expect_equal(gomp$hazard(4), 
                 flexsurv::hgompertz(4, shape = sh, rate = r))     
    
    ## cumhazard
    expect_equal(gomp$cumhazard(4), 
                 flexsurv::Hgompertz(4, shape = sh, rate = r))    
    
    ## random
    set.seed(101)
    r1 <- gomp$random()
    set.seed(101)
    r2 <- flexsurv::rgompertz(1, shape = sh, rate = r)
    expect_equal(r1, r2)
    
    ## Truncated random
    r <- replicate(10, gomp$trandom(2, 9))
    expect_true(all(r >= 2))
    expect_true(all(r <= 9))    
  }
  
  test_gompertz(.05, .5)   # shape > 0
  test_gompertz(0, .5)   # shape = 0
})

# Log-logistic distribution ----------------------------------------------------
test_that("loglogistic", {
  LogLogistic <- module$loglogistic
  sh <- 1; sc <- .5
  llogis <- new(LogLogistic, shape = sh, scale = sc)
  
  # pdf
  expect_equal(llogis$pdf(2), 
              flexsurv::dllogis(2, shape = sh, scale = sc))    
  
  # cdf
  expect_equal(llogis$cdf(2), 
              flexsurv::pllogis(2, shape = sh, scale = sc))      

  # quantile
  expect_equal(llogis$quantile(.34), 
              flexsurv::qllogis(.34, shape = sh, scale = sc))   

  # hazard
  expect_equal(llogis$hazard(6), 
              flexsurv::hllogis(6, shape = sh, scale = sc))     
  
  # cumhazard
  expect_equal(llogis$cumhazard(6), 
              flexsurv::Hllogis(6, shape = sh, scale = sc))   

  # random
  set.seed(101)
  r1 <- llogis$random()
  set.seed(101)
  r2 <- flexsurv::rllogis(1, shape = sh, scale = sc)
  expect_equal(r1, r2)
  
  ## Truncated random
  r <- replicate(10, llogis$trandom(2, 9))
  expect_true(all(r >= 2))
  expect_true(all(r <= 9))   
})

# Generalized gamma distribution -----------------------------------------------
test_that("gengamma", {
  
  GeneralizedGamma <- module$gengamma

  test_gengamma <- function(m, s, q){
    
    gengamma <- new(GeneralizedGamma, mu = m, sigma = s, Q = q)
    
    ## pdf
    expect_equal(gengamma$pdf(2), 
                flexsurv::dgengamma(2, mu = m, sigma = s, Q = q))      
    
    ## cdf
    expect_equal(gengamma$cdf(3), 
                flexsurv::pgengamma(3, mu = m, sigma = s, Q = q))         

    ## quantile
    if (q !=0){
      expect_equal(gengamma$quantile(.3), 
                  flexsurv::qgengamma(.3, mu = m, sigma = s, Q = q))
      expect_equal(gengamma$quantile(.8), 
                  flexsurv::qgengamma(.8, mu = m, sigma = s, Q = q))  
      expect_equal(gengamma$quantile(gengamma$cdf(3)), 3)
    } else{ # flexsurv does not seem to be correct with Q = 0.
      expect_equal(gengamma$quantile(gengamma$cdf(3)), 3)
    }

    ## hazard
    expect_equal(gengamma$hazard(1.8), 
                flexsurv::hgengamma(1.8, mu = m, sigma = s, Q = q))      
    
    ## cumhazard
    expect_equal(gengamma$cumhazard(1.2), 
                flexsurv::Hgengamma(1.2, mu = m, sigma = s, Q = q))
    
    if (q == 0){ # lognormal case
      ## random
      set.seed(101)
      r1 <- gengamma$random()
      set.seed(101)
      r2 <- flexsurv::rgengamma(1, mu = m, sigma = s, Q = q)
      expect_equal(r1, r2)   
    }
  
    ## Truncated random
    r <- replicate(10, gengamma$trandom(2, 9)) 
    expect_true(all(r >= 2))
    expect_true(all(r <= 9)) 
  }
  
  test_gengamma(m = 2, s = 1.5, q = -2) # Q < 0
  test_gengamma(m = 5, s = 2, q = 0) # Q = 0
  test_gengamma(2, 1.1, 2) # Q > 0
})

# Spline survival distribution -------------------------------------------------
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

test_that("survspline", {

  SurvSpline <- module$survspline

  # Scale is log hazard
  test_survspline1 <- function(gamma, knots, timescale){
    
    spline <- new(SurvSpline, gamma = gamma, knots = knots,
                  scale = "log_hazard", timescale = timescale,
                  cumhaz_method = "quad", step = -1, 
                  random_method = "invcdf")
    
    ## hazard
    expect_equal(spline$hazard(2), 
                exp(R_linear_predict(2, gamma, knots, timescale)))  
    
    ## cumhazard
    R_hazard <- function(t) {
      exp(R_linear_predict(t, gamma = gamma, knots = knots,
                          timescale = timescale))
    }
    R_cumhazard <- function(t){
      stats::integrate(R_hazard, 0, t)$value
    }    
    expect_equal(spline$cumhazard(2), 
                 R_cumhazard(2),
                 tolerance = .001, scale = 1)   
    
    ## cdf
    expect_equal(spline$cdf(.8),
                 1 - exp(-spline$cumhazard(.8)))
    
    ## pdf
    R_cdf <- function(t){
      return(1 - exp(-R_cumhazard(t)))
    }
    expect_equal(spline$pdf(0), 
                 0)
    expect_equal(spline$pdf(.5), 
                 (1 - R_cdf(.5)) * R_hazard(.5))
    expect_equal(spline$pdf(.5), 
                  numDeriv::grad(R_cdf, .5))   
    
    ## quantile
    expect_equal(spline$quantile(spline$cdf(.55)), .55, tolerance = .001, scale = 1)
    
    ## random
    expect_error(spline$random(),
                 NA)
    
    ## truncated random
    r <- replicate(10, spline$trandom(2, 2.5)) 
    expect_true(all(r >= 2))
    expect_true(all(r <= 2.5)) 
    
  }

  g = c(-1.2, 1.3, .07)
  k = c(.19, 1.7, 6.7)
  test_survspline1(gamma = g, knots = k, timescale = "identity")
  test_survspline1(gamma = g, knots = k, timescale = "log")
  
  # scale is log cummulative hazard, log odds, or normal
  test_survspline2 <- function(flexsurv_scale, timescale,
                                gamma = NULL, knots = NULL){
    if (is.null(gamma) | is.null(knots)){
      spl_fit <- flexsurvspline(Surv(recyrs, censrec) ~ 1, 
                          data = bc, k = 1, timescale = timescale,
                          scale = flexsurv_scale)
      gamma <- spl_fit$res.t[, "est"]
      knots <- spl_fit$knots
    }
    hesim_scale <- switch(flexsurv_scale,
                          hazard = "log_cumhazard",
                          odds = "log_cumodds",
                          normal = "inv_normal"
                          )
    
    spline <- new(SurvSpline, gamma = gamma, knots = knots,
                  scale = hesim_scale, timescale = timescale,
                  cumhaz_method = "quad", step = -1, random_method = "invcdf")
    
    ### pdf
    expect_equal(flexsurv::dsurvspline(5, gamma = gamma, knots = knots, 
                             scale = flexsurv_scale, timescale = timescale),
                 spline$pdf(5))
    expect_equal(spline$pdf(0),
                 0)
    expect_equal(spline$pdf(-5),
                 0)
    
    ### cdf
    expect_equal(flexsurv::psurvspline(5, gamma = gamma, knots = knots, 
                             scale = flexsurv_scale, timescale = timescale),
                 spline$cdf(5))    
    expect_equal(spline$cdf(0),
                 0)
    expect_equal(spline$cdf(-5),
                 0)   
    
    ### quantile
    expect_equal(flexsurv::qsurvspline(.8, gamma = gamma, knots = knots, 
                             scale = flexsurv_scale, timescale = timescale),
                 spline$quantile(.8),
                  tolerance = .001, scale = 1)      
    
    expect_equal(spline$quantile(-2),
                  NaN)
    expect_equal(spline$quantile(0),
                  -Inf)    
    expect_equal(spline$quantile(1),
                  Inf)        
    
    ### hazard
    expect_equal(flexsurv::hsurvspline(3, gamma = gamma, knots = knots, 
                             scale = flexsurv_scale, timescale = timescale),
                 spline$hazard(3))     
    expect_equal(spline$hazard(0),
                 0)
    expect_equal(spline$hazard(-5),
                 0)    
    
    ### cumhazard
    expect_equal(flexsurv::Hsurvspline(3, gamma = gamma, knots = knots, 
                             scale = flexsurv_scale, timescale = timescale),
                 spline$cumhazard(3))     
    expect_equal(spline$cumhazard(0),
                 0)
    expect_equal(spline$cumhazard(-5),
                 0)   
    
    ### random
    set.seed(12)
    r1 <- rsurvspline(1, gamma = gamma, knots = knots, 
                      scale = flexsurv_scale, timescale = timescale)
    set.seed(12)
    r2 <- spline$random()
    expect_equal(r1, r2, tolerance = .001, scale = 1)
    
    ## truncated random
    r <- replicate(10, spline$trandom(2, 5)) 
    expect_true(all(r >= 2))
    expect_true(all(r <= 5))     
  }
  test_survspline2(flexsurv_scale = "hazard", timescale = "log")
  test_survspline2(flexsurv_scale = "hazard", timescale = "identity")
  test_survspline2(flexsurv_scale = "odds", timescale = "log")
  test_survspline2(flexsurv_scale = "odds", timescale = "identity")
  test_survspline2(flexsurv_scale = "normal", timescale = "log")
  test_survspline2(flexsurv_scale = "normal", timescale = "identity",
                   gamma = c(-1.2, 1.3, .07), knots = c(.19, 1.7, 6.7))
})

# Survival Fractional Polynomials ----------------------------------------------
test_that("FracPoly", {
  
  FracPoly <- module$fracpoly
  
  test_fracpoly <- function(gamma, powers){
    R_hazard <- function(t){
      return(exp(c(gamma %*% t(cbind(1, bfp(t, powers))))))
    }
    R_cumhazard <- function(t){
      stats::integrate(R_hazard, 0, t)$value
    }
    R_linear_predict <- function(t){
      c(gamma %*% t(cbind(1, bfp(t, powers))))
    }    
    
    fp <- new(FracPoly, gamma = gamma, powers = powers, 
              cumhaz_method = "quad", step = -1,
              random_method = "invcdf")
    
    ## hazard
    expect_equal(fp$hazard(4),
                 exp(R_linear_predict(4)))    
    
    ## cumhazard
    expect_equal(fp$cumhazard(2), 
                 R_cumhazard(2),
                 tolerance = .001, scale = 1)
    expect_equal(fp$cumhazard(0), 
                 0) 
    
    ## cdf   
    expect_equal(fp$cdf(2.5), 
                 1 - exp(-R_cumhazard(2.5)),
                 tolerance = .001, scale = 1)  
    
    ## pdf
    R_cdf <- function(t){
      1 - exp(-R_cumhazard(t))
    }  
    expect_equal(fp$pdf(.8), 
                 numDeriv::grad(R_cdf, .8),
                 tolerance = .001, scale = 1)
    
    ## quantile
    expect_equal(fp$quantile(fp$cdf(.45)),
                 .45,
                 tol = .001, scale = 1) 
    
    ## random
    ### Using numeric quantile when max_x_ is infinity
    set.seed(101)
    r1 <- fp$random()
    set.seed(101)
    r2 <- fp$quantile(runif(1, 0, 1))
    expect_equal(r1, r2)
    
    ### Using discrete cumulative hazard otherwise
    fp$max_x_ <- 30
    expect_error(fp$random(), NA) # See manual tests for correctness
    
    ## truncated random
    ### With max_x_ != infinity
    r <- replicate(10, fp$trandom(2, 5)) 
    expect_true(all(r >= 2 & r <= 5))
    
    ### With max_x_ == infinity
    fp$max_x_ <- Inf
    r <- replicate(10, fp$trandom(2, 5)) 
    expect_true(all(r >= 2 & r <= 5))
    
   } # end test_fracpoly
  
  test_fracpoly(gamma = c(-1.2, -.567, 1.15),
                powers = c(1, 0))
  test_fracpoly(gamma = c(.2, .2),
                powers = 0) 
  test_fracpoly(gamma = c(.2, .2, .2),
                powers = c(1, 1))   
  test_fracpoly(gamma = c(.2, .2, .5, .5, .2, .2),
                powers = c(1, 1, 2, 2, 1))
  
  # Equivalence of fractional polynomial and gompertz
  rate <- 1/5; shape <- 1
  fp <- new(FracPoly, gamma = c(log(rate), shape), powers = 1,
            cumhaz_method = "quad", step = -1, random_method = "invcdf")
  expect_equal(fp$hazard(5),
               flexsurv::hgompertz(5, shape = shape, rate = rate))
  expect_equal(fp$cdf(2),
               flexsurv::pgompertz(2, shape = shape, rate = rate))
  
  # Equivalence of fractional polynomial and weibull
  a0 <- -2; a1 <- 1
  fp <- new(FracPoly, gamma = c(a0, a1), powers = 0,
            cumhaz_method = "quad", step = -1, random_method = "invcdf")
  expect_equal(fp$hazard(3), 
               hweibullNMA(3, a0 = a0, a1 = a1))
  expect_equal(fp$cdf(4), 
               pweibullNMA(4, a0 = a0, a1 = a1))
})

# Point mass (fixed) distribution ----------------------------------------------
test_that("fixed", {
  PointMass <- module$point_mass
  est <- 3
  pm <- new(PointMass, est = est)
  
  # pdf
  expect_equal(pm$pdf(est), 1)
  expect_equal(pm$pdf(est - 1), 0)
  
  # cdf
  expect_equal(pm$cdf(est), 1) 
  expect_equal(pm$cdf(est - 1), 0) 
  expect_equal(pm$cdf(est + 1), 1) 
  
  # quantile
  expect_error(pm$quantile(.025))
  
  # hazard
  expect_error(pm$hazard(1))
  
  # cumhazard
  expect_error(pm$cumhazard(4))
  
  # random
  expect_equal(pm$random(), est)
  
  # Truncated random
  expect_equal(pm$trandom(1, Inf), est + 1)
  expect_equal(pm$trandom(1, 2), 2)
  expect_equal(pm$trandom(5, Inf), est + 5)
})

# Truncated normal distribution ------------------------------------------------
test_that("rtruncnorm", {
  n <- 1000
  mu <- 50; sigma <- 10; lower <- 25; upper <- 60
  
  #rtruncnorm from hesim
  set.seed(10)
  samp1 <- replicate(n, hesim:::C_test_rtruncnorm(mu, sigma, lower, upper))
  
  # rtruncnorm from truncnorm package
  set.seed(10)
  samp3 <- truncnorm::rtruncnorm(n, lower, upper, mu, sigma)
  expect_equal(samp1, samp3)
})
