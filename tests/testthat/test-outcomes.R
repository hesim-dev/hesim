context("outcomes.R unit tests")
rm(list = ls())

# Test incr_effect -------------------------------------------------------------
# See unit tests in test-cea.R

# Test surv_quantile -----------------------------------------------------------
test_that("surv_quantile", {
  t <- seq(0, 10, by = .01)
  surv1 <- seq(1, .3, length.out = length(t))
  surv2 <- seq(1, .2, length.out = length(t))
  strategies <- c("Strategy 1", "Strategy 2")
  surv <- data.table(strategy = rep(strategies, each = length(t)),
                     t = rep(t, 2), 
                     surv = c(surv1, surv2))
  quantiles <- surv_quantile(surv, probs = c(.4, .5), t = "t",
                             surv_cols = "surv", by = "strategy")
  row <- which(surv[strategy == "Strategy 1", surv <= 1 - .4])[1]
  expect_true(inherits(quantiles, "data.table"))
  expect_equal(surv[strategy == "Strategy 1"][row, t],
               quantiles[strategy == "Strategy 1" & prob == .4, quantile_surv])
  
  # Check errors
  expect_error(surv_quantile(surv, probs = 1.1, t = "t",
                             surv_cols = "surv", by = "strategy"))
  
  # Check NA handling
  surv2 <- seq(1, .8, length.out = length(t))
  surv <- data.table(strategy = rep(strategies, each = length(t)),
                   t = rep(t, 2), 
                   surv = c(surv1, surv2))
  quantiles <- surv_quantile(surv, probs = c(.4, .5), t = "t",
                             surv_cols = "surv", by = "strategy")
  expect_equal(quantiles$strategy, rep(strategies, 2))
  expect_equal(quantiles$prob, rep(c(.4, .5), each = 2))
  expect_true(all(is.na(quantiles[strategy == "Strategy 2", quantile_surv])))
  
  surv1 <- seq(1, .8, length.out = length(t))
  surv <- data.table(strategy = rep(strategies, each = length(t)),
                   t = rep(t, 2), 
                   surv = c(surv1, surv2))  
  quantiles <- surv_quantile(surv, probs = c(.4, .5), t = "t",
                             surv_cols = "surv", by = "strategy")
  expect_true(all(is.na(quantiles[, quantile_surv])))  
})