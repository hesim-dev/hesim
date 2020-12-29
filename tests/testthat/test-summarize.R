context("outcomes.R unit tests")
rm(list = ls())
library("data.table")

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

# Test summary.ce() ------------------------------------------------------------
# Setup "ce" object
n_grps <- n_dr <- n_sample <- n_strategies <- 2
N <- n_grps * n_dr * n_sample * n_strategies

## Costs
drug_costs <- runif(N, 10000, 40000)
medical_costs <- runif(N, 5000, 10000)
costs <- data.table(
  category = rep(c("drugs", "medical"), each = N),
  dr = rep(rep(c(0, .03), each = n_grps * n_sample * n_strategies),
           2),
  sample = rep(rep(c(1, 2), each = n_grps * n_strategies),
               times = n_dr * n_sample),
  strategy_id = rep(rep(c(1, 2), each = n_grps),
                    times = n_dr * n_sample * n_strategies),
  grp_id = rep(c(1, 2), times = N),
  costs = c(drug_costs, medical_costs)
)
costs_total <- costs[, .(costs = sum(costs)),
                     by = c("dr", "sample", "strategy_id", "grp_id")]
costs_total[, category := "total"]
costs <- rbind(costs, costs_total)

## QALYs
qalys <- costs[category == "drugs"][, c("category", "costs") := NULL]
qalys[, "qalys"] <- runif(N, 3, 5)

## 'ce' object
ce <- list(costs = costs, qalys = qalys)
class(ce) <- "ce"


# Run tests
test_that("summary.ce() returns the correct number or rows", {
  expect_equal(
    nrow(summary(ce, tidy = TRUE)),
    n_grps * n_strategies *  n_dr * 4
  )
  expect_equal(
    nrow(summary(ce, tidy = FALSE)),
    n_grps *  n_dr * 4
  )
})

test_that("summary.ce() correctly passed 'strategy_names' argument", {
  snames <- c("s1", "s2")
  expect_true(
    all(snames %in% colnames(summary(ce, strategy_names = snames)))
  )
})

test_that("summary.ce() correctly passed 'grp_names' argument", {
  gnames <- c("g1", "g2")
  x <- summary(ce, grp_names = gnames)
  expect_equal(as.character(unique(x$grp)), 
               gnames)
})


