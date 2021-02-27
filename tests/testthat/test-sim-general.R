context("sim-general.R unit tests")
library("data.table")

# Test survival() function -----------------------------------------------------
# Predict survival
## Fit models
onc3_pfs_os <- as_pfs_os(onc3, patient_vars = c("patient_id", "female",
                                                "strategy_name"))
fit_pfs <- coxph(Surv(pfs_time, pfs_status) ~ strategy_name + female,
                 data = onc3_pfs_os)
fit_os <- coxph(Surv(os_time, pfs_status) ~ strategy_name + female,
                data = onc3_pfs_os)

## Prediction
newdat <- data.table(
  sample = 1,
  strategy_id = rep(1:3, 2),
  strategy_name = c("SOC", "New 1", "New 2"),
  patient_id = rep(1:2, each = 3),
  female = rep(c(1, 0), each = 3),
  grp_id = 1
)
times <- seq(0, 14, 1/12)
predict_survival <- function(object, newdata, times) {
  surv <- summary(survfit(object, newdata = newdata, se.fit = FALSE),
                  t = times)
  pred <- newdata[rep(seq_len(nrow(newdata)), each = length(times)), ]
  pred[, sample := 1] # Point estimates only in this example
  pred[, time := rep(surv$time, times = nrow(newdata))]
  pred[, survival := c(surv$surv)]
  return(pred[, ])
}
pfs <- predict_survival(fit_pfs, newdata = newdat, times = times)
os <- predict_survival(fit_os, newdata = newdat, times = times)
surv <- rbind(
  as.data.table(pfs)[, curve := 1L],
  as.data.table(os)[, curve := 2L]
)

# Run tests
test_that("$survival() constructs a survival object", {
  s <- survival(surv, t = "time")
  expect_true(inherits(s, "survival"))
})

test_that("$survival() throws an error if number of values within an ID variable is wrong", {
  s <- surv[!(patient_id == 1 & strategy_id == 1)]
  expect_error(
    survival(s, t = "time"),
    paste0("The number of rows in 'data' must be equal to the product of the number ",
           "of unique values of the 'sample', 'strategy_id', 'patient_id' 'grp_id', ",
           "'curve', and 't' columns.")
  )
})