context("Random number generation")
library("msm")

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

  
  # test equal
  expect_equal(samp1, samp2)
})
