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

# State probabilities  ---------------------------------------------------------
sample <- strategy_id <-  rep(1, 5) - 1
patient_id <- c(rep(1, 3), rep(2, 2)) - 1
from <- c(1, 2, 1, 1, 2) - 1
to <- c(2, 1, 3, 2, 3) - 1
final <- c(0, 0, 1, 0, 1)
time_start <- c(0, 2, 5, 0, 3)
time_stop <- c(2, 5, 7, 3, 8)

disprog <- data.frame(sample = sample, 
                      strategy_id = strategy_id,
                      patient_id = patient_id,
                      from = from,
                      to = to,
                      final = final,
                      time_start = time_start, 
                      time_stop = time_stop)

test_that("C_ctstm_indiv_stateprobs", {
  # Time 2.5
  stprobs <- data.table(hesim:::C_ctstm_indiv_stateprobs(disprog, t = 2.5, n_samples = 1,
                                              n_strategies = 1, n_states = 3,
                                              n_patients = 2))
  expect_equal(stprobs[state_id == 1, prob], .5)
  expect_equal(stprobs[state_id == 2, prob], 0)
  
  # Time 3
  stprobs <- data.table(hesim:::C_ctstm_indiv_stateprobs(disprog, t = 3, n_samples = 1,
                                              n_strategies = 1, n_states = 3,
                                              n_patients = 2))
  expect_equal(stprobs[state_id == 1, prob], 1)
  
  # Time 10
  stprobs <- data.table(hesim:::C_ctstm_indiv_stateprobs(disprog, t = 10, n_samples = 1,
                                              n_strategies = 1, n_states = 3,
                                              n_patients = 2))
  expect_equal(stprobs[state_id == 1, prob], 0)
  expect_equal(stprobs[state_id == 2, prob], 1)
})

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

test_that("CtstmTrans", {
  # hazard
  hesim_hazard <- transmod$hazard(3)
  expect_equal(hesim_hazard[trans == 1][1]$hazard,
               summary(fits[[1]], type = "hazard", t = 3, ci = FALSE)[[1]][1, "est"])
  expect_equal(hesim_hazard[trans == 2][1]$hazard,
               summary(fits[[2]], type = "hazard", t = 3, ci = FALSE)[[1]][1, "est"])

  # cumulative hazard
  hesim_cumhazard <- transmod$cumhazard(5)
  expect_equal(hesim_cumhazard[trans == 1][1]$cumhazard,
               summary(fits[[1]], type = "cumhaz", t = 5)[[1]][1, "est"])
  expect_equal(hesim_cumhazard[trans == 2][1]$cumhazard,
               summary(fits[[2]], type = "cumhaz", t = 5)[[1]][1, "est"])
  
  # State probabilities
  expect_error(transmod$sim_stateprobs(t = c(0, 1, 2, 3)),
               NA)
  
  # Errors
  expect_error(create_CtstmTrans(fits, data = fits_data, trans_mat = tmat,
                                start_ages = rep(1, nrow(dt_patients)),
                                point_estimate = TRUE),
               NA)
  expect_error(create_CtstmTrans(fits, data = fits_data, trans_mat = tmat,
                                start_ages = rep(1, nrow(dt_patients)),
                                death_state = nrow(tmat),
                                point_estimate = TRUE),
               NA)
  expect_error(create_CtstmTrans(fits, data = fits_data, trans_mat = tmat,
                                start_ages = rep(1, nrow(dt_patients) - 1),
                                point_estimate = TRUE))
  expect_error(create_CtstmTrans(fits, data = fits_data, trans_mat = tmat,
                                death_state = nrow(tmat) + 1,
                                point_estimate = TRUE))
})  

test_that("IndivCtstm", {
  ictstm <- IndivCtstm$new(trans_model = transmod)
  
  # Errors 
  expect_error(ictstm$sim_stateprobs())
  transmod2 <- transmod$clone()
  transmod2$trans_mat <- matrix(1)
  expect_error(transmod2$sim_stateprobs(t = 2))
  transmod2$trans_mat <-1
  expect_error(transmod2$sim_stateprobs(t = 2))
  transmod2$trans_mat <- matrix(seq(1, 6), nrow = 2)
  expect_error(transmod2$sim_stateprobs(t = 2))
  
  # Base case simulation
  ## Simulate disease progression
  set.seed(101)
  disprog <- ictstm$sim_disease()$disease_prog_
  
  ## Time from first state
  set.seed(101)
  time_1a <- rexp(1, rate = fits[[1]]$res[, "est"])
  time_1b <- rexp(1, rate = fits[[2]]$res[, "est"])
  time1 <- min(time_1a, time_1b)
  state1 <- which.min(c(time_1a, time_1b)) + 1
  expect_equal(disprog[1, time_stop], time1)
  expect_equal(disprog[1, to], state1)
  
  ### Time from second state
  time2 <- rexp(1, rate = fits[[3]]$res[, "est"])
  expect_equal(disprog[2, time_stop],  time1 + time2)
  
  ## Simulate state probabilities
  stprobs <- ictstm$sim_stateprobs(t = c(0, 1, 2, 3))$stateprobs_
  
  expect_true(all(stprobs[t == 0 & state_id == 1, prob] ==  1))
  expect_true(all(stprobs[t == 0 & state_id != 1, prob] ==  0))
  
  # Simulation other cases 
  ## Maximum time = 2
  disprog <- ictstm$sim_disease(max_t = 2)$disease_prog_
  expect_equal(ictstm$stateprobs_, NULL)
  expect_equal(max(disprog$time_stop), 2)
  expect_true(all(disprog[final == 1 & time_stop == 2, to] == 1)) # All should have remained in initial state
  
  ## Maximum age = 43 (i.e., max_t = 5)
  disprog <- ictstm$sim_disease(max_age = 43)$disease_prog_
  expect_true(all(disprog[time_stop == 5 & final == 1, to] == 3)) # All should have moved to death state
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
fit <- flexsurvreg(Surv(time/365.25, status) ~ factor(trans), data = msebmt,
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

test_that("CtstmTrans", {
  # hazard
  hesim_hazard <- transmod$hazard(3)
  expect_equal(hesim_hazard[trans == 1][1]$hazard,
               summary(fit, type = "hazard", t = 3, ci = FALSE)[[1]][1, "est"])
  expect_equal(hesim_hazard[trans == 4][1]$hazard,
               summary(fit, type = "hazard", t = 3, ci = FALSE)[[4]][1, "est"])

  # Cumulative hazard
  hesim_cumhazard <- transmod$cumhazard(3)
  expect_equal(hesim_cumhazard[trans == 1][1]$cumhazard,
               summary(fit, type = "cumhaz", t = 3, ci = FALSE)[[1]][1, "est"])
  expect_equal(hesim_cumhazard[trans == 4][1]$cumhazard,
               summary(fit, type = "cumhaz", t = 3, ci = FALSE)[[4]][1, "est"])
})

test_that("IndivCtstm", {
  ictstm <- IndivCtstm$new(trans_model = transmod)
  
  # Simulate disease progression
  expect_error(ictstm$sim_disease()$disease_prog_, NA)
  disprog <- ictstm$sim_disease()$disease_prog_
  expect_true(is.data.table(disprog))
  
  # Simulate state probabilities
  stprobs <- ictstm$sim_stateprobs(t = c(0, 1, 2, 3))$stateprobs_
  expect_true(all(stprobs[t == 0 & state_id == 1, prob] ==  1))
  expect_true(all(stprobs[t == 0 & state_id != 1, prob] ==  0))
})