context("Random number generation")
library("msm")
library("truncnorm")

# Test rcat -------------------------------------------------------------------
test_that("rcat", {

  # does sampling work for a vector?
  p <- c(.2, .3, .5)
  samp <- rcat(10000, p)
  expect_equal(as.numeric(prop.table(table(samp))), p, tolerance = .03)
  
  # does sampling work for a matrix with rows less than rows in n?
  p1 <- c(.3, .5, .2)
  p2 <- c(.1, .1, .8)
  pmat <- matrix(c(p1, p2), nrow = 2, byrow = TRUE)
  samp <- rcat(10000, pmat)
  samp1 <- samp[seq(1, length(samp), by = 2)]
  samp2 <- samp[seq(2, length(samp), by = 2)]
  expect_equal(as.numeric(prop.table(table(samp1))), p1, tolerance = .03)
  expect_equal(as.numeric(prop.table(table(samp2))), p2, tolerance = .03)
  
  # should be equal to rmultinom  
  p <- c(.2, .5, .3, .5, .1, .4)
  n <- 1000
  pmat <- matrix(rep(p, n/2), nrow = n, ncol = 3, byrow = TRUE)
  
  ## rcat
  set.seed(100)
  samp1 <- rcat(n, pmat)
  prop.table(table(samp1))
  
  ## rmultinom
  set.seed(100)
  samp2 <- t(apply(pmat, 1, rmultinom, n = 1, size = 1))
  samp2 <- apply(samp2, 1, function(x) which(x == 1))
  prop.table(table(samp2))
  
  ## test equal
  expect_equal(samp1, samp2)
  
  # n must be positive
  expect_error(rcat(n = -1, prob = c(.2, .8)))
})

# Test rpwexp -----------------------------------------------------------------
test_that("rpwexp", {
  
  # does sampling work for a vector?
  n <- 10000
  r <- c(.6, 1.2, 1.3)
  t <- c(0, 10, 15)
  samp1 <- rpwexp(n, rate = r, time = t)
  samp2 <- msm::rpexp(n, r, t)
  expect_equal(mean(samp1), mean(samp2), tol = .1, scale = 1)
  
  # does sampling work for a matrix with rows less than rows in n?
  n <- 20000
  r1 <- c(.6, 1.2, 1.3)
  r2 <- c(1.3, 1.8, 2.0)
  rmat <- matrix(c(r1, r2), nrow = 2, byrow = TRUE)
  samp1 <- rpwexp(n, rate = rmat, time = t)
  samp2.1 <- msm::rpexp(n, rate = r1, t = t)
  samp2.2 <- msm::rpexp(n, rate = r2, t = t)
  expect_equal(mean(samp1[seq(1, length(samp1), by = 2)]), mean(samp2.1), tol = .1, scale = 1)
  
  # does sampling work for a matrix with rows equals to n?
  n <- 10000
  rmat <- matrix(c(.6, 1.2, 1.3), ncol = 3, nrow = n, byrow = TRUE)
  samp1 <- rpwexp(n, rate = rmat, time = t)
  samp2 <- msm::rpexp(n, rate = rmat[1, ], t = t)
  expect_equal(median(samp1), median(samp2), tol = .1, scale = 1)
  
  # check errors
  expect_error(rpwexp(2, rate = c(.2, 1.2, 2.0), time = c(.2, 1)))
  expect_error(rpwexp(2, rate = c(.2, 1.2, 2.0), 
                      time = matrix(c(.2, 1), ncol = 2, byrow = TRUE)))
  expect_error(rpwexp(n = -2))
  expect_error(rpwexp(n = 2, rmat))
})

# Test rdirichlet_mat ---------------------------------------------------------
rdirichlet <- function(n, alpha){
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}

# check accuracy
n <- 10000
alpha <- matrix(c(100, 200, 500, 50, 70, 75), ncol = 3, nrow = 2, byrow = TRUE)
samp1 <- rdirichlet_mat(n, alpha)
samp2.1 <- rdirichlet(n, alpha[1, ])
mean1 <- apply(samp1, c(1, 2), mean)
mean2.1 <- apply(samp2.1, 2, mean)
expect_equal(mean1[1, ], mean2.1, tolerance = .02, scale = 1)

# works with vector
expect_error(rdirichlet_mat(n = 1, alpha[1, ]), NA)

# check errors
expect_error(rdirichlet_mat(n = -1, alpha))
expect_error(rdirichlet_mat(n = 1, array(alpha, dim = c(2, 3, 1))))

# Test c++ function rtruncnormC -----------------------------------------------
test_that("rtruncnormC", {
  n <- 1000
  mu <- 50; sigma <- 10; lower <- 25; upper <- 60
  
  #rtruncnormC from hesim
  set.seed(10)
  samp1 <- replicate(n, hesim:::rtruncnormC(mu, sigma, lower, upper))
  
  # rtnorm from msm package
  set.seed(10)
  samp2 <- rtnorm(n, mu, sigma, lower, upper)
  
  # rtruncnorm from truncnorm package
  set.seed(10)
  samp3 <- truncnorm::rtruncnorm(n, lower, upper, mu, sigma)
  expect_equal(samp1, samp3)
})

