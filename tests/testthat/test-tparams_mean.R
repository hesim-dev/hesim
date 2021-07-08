context("test-tparqms_mean.R unit tests")

hesim_dat <- hesim_data(
  strategies = data.frame(strategy_id = c(1, 2)),
  patients = data.frame(patient_id = c(1, 2)),
  states = data.frame(
    state_id = c(1, 2, 3),
    state_name = c("state1", "state2", "state3")
  )
)

# Time-invariant cost model
cost_tbl <- stateval_tbl(
  data.frame(strategy_id = hesim_dat$strategies$strategy_id,
             mean = c(5000, 3000),
             se = c(200, 100)
  ),
  dist = "gamma"
)
costmod <- create_StateVals(cost_tbl, n = 2, hesim_data = hesim_dat)

# Time-varying cost model
cost_tbl_tv <- stateval_tbl(
  data.frame(strategy_id = hesim_dat$strategies$strategy_id,
             time_start = c(0, 0, 5, 5),
             mean = c(5000, 3000, 5000 * 2, 3000 * 2),
             se = c(200, 100, 200 * 2, 100 * 2)
  ),
  dist = "gamma"
)
costmod_tv <- create_StateVals(cost_tbl_tv, n = 2, hesim_data = hesim_dat)

# Summary method ---------------------------------------------------------------
test_that("summary.tparams_mean works as expected", {
  ps <- summary(costmod$params, probs = c(.2, .5, .7))
  expect_true(inherits(ps, "data.table"))
  expect_equal(
    colnames(ps),
    c("strategy_id", "patient_id", "state_id",
      "mean", "sd", "20%", "50%", "70%")
  )
  expect_equal(ps$mean, apply(costmod$params$value, 1, mean))
})

test_that("summary.tparams_mean works with prob of lenngth 1", {
  ps <- summary(costmod$params, probs = .5)
  expect_equal(
    colnames(ps),
    c("strategy_id", "patient_id", "state_id",
      "mean", "sd", "50%")
  )
})

test_that("summary.tparams_mean works withtime-varying means", {
  ps <- summary(costmod_tv$params)
  expect_equal(
    colnames(ps),
    c("strategy_id", "patient_id", "state_id",
      "time_id", "time_start", "time_stop",
      "mean", "sd", "2.5%", "97.5%")
  )
})

# Print method -----------------------------------------------------------------
test_that("print.params_surv() works as expected", {
  p <- costmod$params
  expect_output(print(p), "A \"tparams_mean\" object")
  expect_output(print(p), "Summary of means:")
})