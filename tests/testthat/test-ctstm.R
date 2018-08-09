context("ctstm.R unit tests")
library("flexsurv")
library("mstate")
library("data.table")
rm(list = ls())

# Simulation data
dt_strategies <- data.table(strategy_id = c(1, 2, 3))
dt_patients <- data.table(patient_id = seq(1, 3),
                          age = c(45, 50, 60),
                          female = c(0, 0, 1))

# List of models  --------------------------------------------------------------
# Multi-state model
fits <- vector(length = 3, mode = "list")
dat <- data.table(bosms3)
for (i in 1:length(fits)){
  fits[[i]] <- flexsurvreg(Surv(years, status) ~ 1, data = dat[trans == i],
                           dist = "exp")
}
fits <- flexsurvreg_list(fits)
tmat <- rbind(c(NA, 1, 2),
              c(NA, NA, 3),
              c(NA, NA, NA))

# Simulation model
hesim_dat <- hesim_data(strategies = dt_strategies,
                        patients = dt_patients)
fits_data <- expand_hesim_data(hesim_dat)
transmod <- create_CtstmTrans(fits, data = fits_data, trans_mat = tmat,
                              point_estimate = TRUE)

test_that("transmod$hazard", {
  hesim_hazard <- transmod$hazard(3)
  expect_equal(hesim_hazard[trans == 1][1]$hazard,
               summary(fits[[1]], type = "hazard", t = 3, ci = FALSE)[[1]][1, "est"])
  expect_equal(hesim_hazard[trans == 2][1]$hazard,
               summary(fits[[2]], type = "hazard", t = 3, ci = FALSE)[[1]][1, "est"])
})

test_that("transmod$cumhazard", {
  hesim_cumhazard <- transmod$cumhazard(5)
  expect_equal(hesim_cumhazard[trans == 1][1]$cumhazard,
               summary(fits[[1]], type = "cumhaz", t = 5)[[1]][1, "est"])
  expect_equal(hesim_cumhazard[trans == 2][1]$cumhazard,
               summary(fits[[2]], type = "cumhaz", t = 5)[[1]][1, "est"])
})

# Joint model  -----------------------------------------------------------------
# Multi-state model
## Data
data("ebmt4")
tmat <- transMat(x = list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6), c(5, 6),
                c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))
msebmt <- msprep(data = ebmt4, trans = tmat, time = c(NA, "rec", "ae",
                 "recae", "rel", "srv"), status = c(NA, "rec.s", "ae.s", "recae.s",
                  "rel.s", "srv.s"), keep = c("match", "proph", "year", "agecl"))
covs <- c("match", "proph", "year", "agecl")
msebmt <- expand.covs(msebmt, covs, longnames = FALSE)

## Fits
fit <- flexsurvreg(Surv(time, status) ~ factor(trans), data = msebmt,
                           dist = "weibull")

# Simulation
dt_transitions <- create_trans_dt(tmat)
dt_transitions[, trans := transition_id]
hesim_dat <- hesim_data(strategies = dt_strategies,
                        patients = dt_patients,
                        transitions = dt_transitions)
fit_data <- expand_hesim_data(hesim_dat, by = c("strategies", "patients", "transitions"))
transmod <- create_CtstmTrans(fit, data = fit_data, trans_mat = tmat,
                              point_estimate = TRUE)

test_that("transmod$hazard", {
  hesim_hazard <- transmod$hazard(3)
  expect_equal(hesim_hazard[trans == 1][1]$hazard,
               summary(fit, type = "hazard", t = 3, ci = FALSE)[[1]][1, "est"])
  expect_equal(hesim_hazard[trans == 4][1]$hazard,
               summary(fit, type = "hazard", t = 3, ci = FALSE)[[4]][1, "est"])
})