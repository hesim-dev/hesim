context("Random number generation")

# Test rcat -------------------------------------------------------------------
test_that("rcat", {
  p <- c(.2, .5, .3)
  pmat <- matrix(rep(p, 10000), nrow = 10000, ncol = length(p), byrow = TRUE)
  
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
