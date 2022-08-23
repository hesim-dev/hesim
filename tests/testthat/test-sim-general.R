context("sim-general.R unit tests")
library("data.table")

# Test survival() function -----------------------------------------------------
# Predict survival
## Fit models
onc3_pfs_os <- as_pfs_os(onc3, patient_vars = c(
  "patient_id", "female",
  "strategy_name"
))
fit_pfs <- coxph(Surv(pfs_time, pfs_status) ~ strategy_name + female,
  data = onc3_pfs_os
)
fit_os <- coxph(Surv(os_time, pfs_status) ~ strategy_name + female,
  data = onc3_pfs_os
)

## Prediction
newdat <- data.table(
  sample = 1,
  strategy_id = rep(1:3, 2),
  strategy_name = c("SOC", "New 1", "New 2"),
  patient_id = rep(1:2, each = 3),
  female = rep(c(1, 0), each = 3),
  grp_id = 1
)
times <- seq(0, 14, 1 / 12)
predict_survival <- function(object, newdata, times) {
  surv <- summary(survfit(object, newdata = newdata, se.fit = FALSE),
    t = times
  )
  pred <- newdata[rep(seq_len(nrow(newdata)), each = length(times)), ]
  pred[, sample := 1] # Point estimates only in this example
  pred[, time := rep(surv$time, times = nrow(newdata))]
  pred[, survival := c(surv$surv)]
  return(pred[, ])
}
pfs <- predict_survival(fit_pfs, newdata = newdat, times = times)
os <- predict_survival(fit_os, newdata = newdat, times = times)
surv_dt <- rbind(
  as.data.table(pfs)[, curve := 1L],
  as.data.table(os)[, curve := 2L]
)

# Run tests
test_that("$survival() constructs a survival object", {
  s <- survival(surv_dt, t = "time")
  expect_true(inherits(s, "survival"))
})

test_that("$survival() throws an error if number of values within an ID variable is wrong", {
  s <- surv_dt[!(patient_id == 1 & strategy_id == 1)]
  expect_error(
    survival(s, t = "time"),
    paste0(
      "The number of rows in 'data' must be equal to the product of the number ",
      "of unique values of the 'sample', 'strategy_id', 'patient_id' 'grp_id', ",
      "'curve', and 't' columns."
    )
  )
})

# Test sim_stateprobs.survival() -----------------------------------------------
surv <- survival(surv_dt, t = "time")
surv[, curve_name := paste0("curve", curve)]
survw <- dcast(surv,
  sample + strategy_id + patient_id + grp_id + t ~ curve_name,
  value.var = "survival"
)
stprobs <- sim_stateprobs(surv)[, state_name := paste0("state", state_id)]
stprobsw <- dcast(stprobs,
  sample + strategy_id + patient_id + grp_id + t ~ state_name,
  value.var = "prob"
)

test_that("The first health state in sim_stateprobs.survival() has the correct probability", {
  expect_equal(stprobsw$state1, survw$curve1)
})

test_that("The middle health states in sim_stateprobs.survival() have the correct probability", {
  expect_equal(stprobsw$state2, survw$curve2 - survw$curve1)
})

test_that("The final health states in sim_stateprobs.survival() has the correct probability", {
  expect_equal(stprobsw$state3, 1 - survw$curve2)
})

surv2 <- survival(
  data.table(
    sample = 1,
    strategy_id = 1,
    patient_id = 1,
    grp_id = 1,
    curve = rep(c(1, 2), each = 2),
    t = rep(c(.6, .8), 2),
    survival = c(.9, .7, .95, .6)
  )
)

test_that("sim_stateprobs.survival() produces expected warning when curves cross", {
  expect_warning(
    sim_stateprobs(surv2),
    "The survival curves were crossed 1/6 (16.7%) of the time.",
    fixed = TRUE
  )
})

test_that("sim_stateprobs.survival() sets probabilities to zero when curves cross", {
  # Curves cross at one time point
  p <- suppressWarnings(sim_stateprobs(surv2))
  expect_true(p[state_id == 2 & t == 0.8]$prob == 0)

  # Curves cross at all time points
  s2 <- copy(surv2)
  s2[t == .6, survival := 1]
  p <- suppressWarnings(sim_stateprobs(s2))
  expect_true(all(p[state_id == 2]$prob == 0))
})

test_that("sim_stateprobs.survival() ensures probabilities sum to 1 (v1)", {
  s <- copy(surv)
  s[, survival := ifelse(curve == 1, 1, survival)]
  p <- suppressWarnings(sim_stateprobs(s))

  # Probabilities are 0 for states 2 and 3
  expect_true(all(p[state_id > 1]$prob == 0))

  # Probabilities sum to 1
  p_sum <- p[, .(prob = sum(prob)), by = c(
    "sample", "strategy_id", "patient_id",
    "grp_id", "t"
  )]
  expect_true(all(p_sum$prob == 1))
})

test_that("sim_stateprobs.survival() ensures probabilities sum to 1 (v2)", {
  p_sum <- stprobs[, .(prob = sum(prob)), by = c(
    "sample", "strategy_id", "patient_id",
    "grp_id", "t"
  )]
  expect_true(all.equal(p_sum$prob, rep(1, nrow(p_sum))))
})

test_that("sim_stateprobs.survival() sets probabilities in successive states to zero when multiple curves cross", {
  s <- rbind(
    surv2,
    data.table(
      sample = 1, strategy_id = 1, patient_id = 1, grp_id = 1,
      curve = 3, t = c(.6, .8), survival = c(.85, .5)
    )
  )
  p <- suppressWarnings(sim_stateprobs(survival(s)))

  # Multiple curves cross at time 0.8
  expect_true(all(p[state_id %in% c(2, 3) & t == .8]$prob == 0))

  # Only 3rd and 2nd curves cross at time 0.6
  expect_true(all(p[state_id == 3 & t == .6]$prob == 0))

  # Probabilities are correct in final state
  expect_equal(p[state_id == 4 & t == .6]$prob, .05)
  expect_equal(p[state_id == 4 & t == .8]$prob, .3)
})
