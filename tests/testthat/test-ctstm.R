context("ctstm.R unit tests")
library("flexsurv")
library("mstate")
library("data.table")
rm(list = ls())

# Helper functions  ------------------------------------------------------------
pv <- function(z, r, t1, t2) {
  z * ((exp(-r * t1) - exp(-r * t2))/r)
}

# Statistical models to estimate parameters  -----------------------------------
n_samples <- 2

# Utility
beta_params <- mom_beta(ctstm3_exdata$utility$mean, ctstm3_exdata$utility$se)
utility1 <- rbeta(n_samples, shape1 = beta_params$shape1[1], 
                  shape2 = beta_params$shape2[1])
utility2 <- rbeta(n_samples, shape1 = beta_params$shape1[2], 
                  shape2 = beta_params$shape2[2])
utility_params <- params_lm(coefs  = cbind(utility1, utility2))

# Costs
medcosts_params <- params_lm(coefs = cbind(runif(n_samples, 5000, 10000),
                                           runif(n_samples, 4000, 8000)))
drugcosts_mat <- matrix(rep(ctstm3_exdata$costs$drugs$costs, n_samples), 
                         nrow = n_samples, byrow = TRUE)
drugcosts_params <- params_lm(coefs  = drugcosts_mat)
 

# Clock-reset multi-state model
## Separate 
msfit_list <- vector(length = 3, mode = "list")
dat <- data.table(bosms3)
for (i in 1:length(msfit_list)){
  msfit_list[[i]] <- flexsurvreg(Surv(years, status) ~ 1, data = dat[trans == i],
                           dist = "exp")
}
msfit_list <- flexsurvreg_list(msfit_list)

## Joint
data("ebmt4")
tmat_ebmt4 <- transMat(x = list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6), c(5, 6),
                c(), c()), names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death"))
msebmt <- msprep(data = ebmt4, trans = tmat_ebmt4, time = c(NA, "rec", "ae",
                 "recae", "rel", "srv"), status = c(NA, "rec.s", "ae.s", "recae.s",
                  "rel.s", "srv.s"), keep = c("match", "proph", "year", "agecl"))
covs <- c("match", "proph", "year", "agecl")
msebmt <- expand.covs(msebmt, covs, longnames = FALSE)
msfit <- flexsurvreg(Surv(time/365.25, status) ~ factor(trans), data = msebmt,
                           dist = "weibull")

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

# Simulate decision model(s) ---------------------------------------------------
# Construct a model structure 
## Simulation data
dt_strategies <- data.table(strategy_id = c(1, 2, 3))
dt_patients <- data.table(patient_id = seq(1, 3),
                          age = c(45, 50, 60),
                          female = c(0, 0, 1))
dt_states <- data.table(state_id = c(1, 2))
hesim_dat <- hesim_data(strategies = dt_strategies,
                        patients = dt_patients,
                        states = dt_states)

## Cost and utility models
statevals_edata <- expand_hesim_data(hesim_dat, by = c("strategies", "patients",
                                                       "states"))
input_dat <- create_input_data(formula_list(mu = ~ -1 + factor(state_id)), 
                               data = statevals_edata)
utilmod <- StateVals$new(data = input_dat, params = utility_params)
medcostsmod <- StateVals$new(data = input_dat, params = medcosts_params)  
drugcostsmod <- StateVals$new(data = input_dat, params = drugcosts_params)  

## The health state transitions
### With transition specific survival models
msfit_list_data <- expand_hesim_data(hesim_dat)
tmat <- rbind(c(NA, 1, 2),
              c(NA, NA, 3),
              c(NA, NA, NA))
mstate_list <- create_CtstmTrans(msfit_list, data = msfit_list_data, trans_mat = tmat,
                                      point_estimate = TRUE)

test_that("CtstmTrans", {
  # hazard
  hesim_hazard <- mstate_list$hazard(3)
  expect_equal(hesim_hazard[trans == 1][1]$hazard,
               summary(msfit_list[[1]], type = "hazard", t = 3, ci = FALSE)[[1]][1, "est"])
  expect_equal(hesim_hazard[trans == 2][1]$hazard,
               summary(msfit_list[[2]], type = "hazard", t = 3, ci = FALSE)[[1]][1, "est"])

  # cumulative hazard
  hesim_cumhazard <- mstate_list$cumhazard(5)
  expect_equal(hesim_cumhazard[trans == 1][1]$cumhazard,
               summary(msfit_list[[1]], type = "cumhaz", t = 5)[[1]][1, "est"])
  expect_equal(hesim_cumhazard[trans == 2][1]$cumhazard,
               summary(msfit_list[[2]], type = "cumhaz", t = 5)[[1]][1, "est"])
  
  # State probabilities
  expect_error(mstate_list$sim_stateprobs(t = c(0, 1, 2, 3)),
               NA)
  
  # Errors
  expect_error(create_CtstmTrans(msfit_list, data = msfit_list_data, trans_mat = tmat,
                                start_ages = rep(1, nrow(dt_patients)),
                                point_estimate = TRUE),
               NA)
  expect_error(create_CtstmTrans(msfit_list, data = msfit_list_data, trans_mat = tmat,
                                start_ages = rep(1, nrow(dt_patients)),
                                death_state = nrow(tmat),
                                point_estimate = TRUE),
               NA)
  expect_error(create_CtstmTrans(msfit_list, data = msfit_list_data, trans_mat = tmat,
                                start_ages = rep(1, nrow(dt_patients) - 1),
                                point_estimate = TRUE))
  expect_error(create_CtstmTrans(msfit_list, data = msfit_list_data, trans_mat = tmat,
                                death_state = nrow(tmat) + 1,
                                point_estimate = TRUE))
})  

### With a joint model
dt_transitions <- create_trans_dt(tmat_ebmt4)
dt_transitions[, trans := transition_id]
hesim_dat$transitions <- dt_transitions
msfit_data <- expand_hesim_data(hesim_dat, by = c("strategies", "patients", "transitions"))
mstate <- create_CtstmTrans(msfit, data = msfit_data, trans_mat = tmat_ebmt4,
                            point_estimate = TRUE)

test_that("CtstmTrans", {
  # hazard
  hesim_hazard <- mstate$hazard(3)
  expect_equal(hesim_hazard[trans == 1][1]$hazard,
               summary(msfit, type = "hazard", t = 3, ci = FALSE)[[1]][1, "est"])
  expect_equal(hesim_hazard[trans == 4][1]$hazard,
               summary(msfit, type = "hazard", t = 3, ci = FALSE)[[4]][1, "est"])

  # Cumulative hazard
  hesim_cumhazard <- mstate$cumhazard(3)
  expect_equal(hesim_cumhazard[trans == 1][1]$cumhazard,
               summary(msfit, type = "cumhaz", t = 3, ci = FALSE)[[1]][1, "est"])
  expect_equal(hesim_cumhazard[trans == 4][1]$cumhazard,
               summary(msfit, type = "cumhaz", t = 3, ci = FALSE)[[4]][1, "est"])
})

# Simulate outcomes
## With transition specific survival models
test_that("Simulate disease and state probabilities", {
  ictstm <- IndivCtstm$new(trans_model = mstate_list,
                           utility_model = utilmod)
  
  # Errors 
  expect_error(ictstm$sim_stateprobs())
  mstate_list2 <- mstate_list$clone()
  mstate_list2$trans_mat <- matrix(1)
  expect_error(mstate_list2$sim_stateprobs(t = 2))
  mstate_list2$trans_mat <-1
  expect_error(mstate_list2$sim_stateprobs(t = 2))
  mstate_list2$trans_mat <- matrix(seq(1, 6), nrow = 2)
  expect_error(mstate_list2$sim_stateprobs(t = 2))
  
  # Base case simulation
  ## Simulate disease progression
  set.seed(101)
  disprog <- ictstm$sim_disease()$disease_prog_
  
  ## Time from first state
  set.seed(101)
  time_1a <- rexp(1, rate = msfit_list[[1]]$res[, "est"])
  time_1b <- rexp(1, rate = msfit_list[[2]]$res[, "est"])
  time1 <- min(time_1a, time_1b)
  state1 <- which.min(c(time_1a, time_1b)) + 1
  expect_equal(disprog[1, time_stop], time1)
  expect_equal(disprog[1, to], state1)
  
  ### Time from second state
  time2 <- rexp(1, rate = msfit_list[[3]]$res[, "est"])
  expect_equal(disprog[2, time_stop],  time1 + time2)
  
  ## Simulate state probabilities
  stprobs <- ictstm$sim_stateprobs(t = c(0, 1, 2, 3))$stateprobs_
  expect_true(all(stprobs[t == 0 & state_id == 1, prob] ==  1))
  expect_true(all(stprobs[t == 0 & state_id != 1, prob] ==  0))
  
  # Simulation with other cases 
  ## Maximum time = 2
  disprog <- ictstm$sim_disease(max_t = 2)$disease_prog_
  expect_equal(ictstm$stateprobs_, NULL)
  expect_equal(max(disprog$time_stop), 2)
  expect_true(all(disprog[final == 1 & time_stop == 2, to] == 1)) # All should have remained in initial state
  
  ## Maximum age = 43 (i.e., max_t = 5)
  disprog <- ictstm$sim_disease(max_age = 43)$disease_prog_
  expect_true(all(disprog[time_stop == 5 & final == 1, to] == 3)) # All should have moved to death state
})

mstate_list <- create_CtstmTrans(msfit_list, data = msfit_list_data, trans_mat = tmat,
                                 n = n_samples)
test_that("Simulate costs and QALYs", {
  ictstm <- IndivCtstm$new(trans_model = mstate_list,
                           utility_model = utilmod,
                           cost_models = list(medical = medcostsmod, 
                                              drugs = drugcostsmod))
  ictstm$sim_disease()
  
  # Simulate QALYs
  ## dr = .03
  ictstm$sim_qalys(dr = .03)$qalys_
  
  disprog1 <- ictstm$disease_prog_[sample == 1 & strategy_id == 1 & patient_id == 2]
  qalys1 <- ictstm$qalys_[sample == 1 & strategy_id == 1 & patient_id == 2]
  utilvals <- utility_params$coefs[1, disprog1$from] 
  qalys_expected <- sum(pv(utilvals, .03, disprog1$time_start, disprog1$time_stop))
  expect_equal(qalys1$qalys, qalys_expected)
  
  ## dr = 0
  ictstm$utility_model$params$coefs <- matrix(1, nrow = n_samples, ncol = nrow(dt_states))
  qalys <- ictstm$sim_qalys(dr = 0)$qalys_
  expect_equal(ictstm$disease_prog_[final == 1][3, time_stop],
               qalys[3, qalys])
  
  # Simulate costs
   costs <- ictstm$sim_costs(dr = c(0, .03))$costs_
   expect_equal(unique(costs$category), c("medical", "drugs"))
   expect_equal(unique(costs$dr), c(0, .03))
   
  # Summarize costs and QALYs
  ce_summary <- ictstm$summarize()
  expect_true(inherits(ce_summary, "ce"))
  ictstm$sim_disease()
  expect_error(ictstm$summarize())
  ictstm$sim_costs()
  expect_error(ictstm$summarize())
})


## With a joint survival model
test_that("IndivCtstm", {
  ictstm <- IndivCtstm$new(trans_model = mstate)
  
  # Simulate disease progression
  expect_error(ictstm$sim_disease()$disease_prog_, NA)
  disprog <- ictstm$sim_disease()$disease_prog_
  expect_true(is.data.table(disprog))
  
  # Simulate state probabilities
  stprobs <- ictstm$sim_stateprobs(t = c(0, 1, 2, 3))$stateprobs_
  expect_true(all(stprobs[t == 0 & state_id == 1, prob] ==  1))
  expect_true(all(stprobs[t == 0 & state_id != 1, prob] ==  0))
})


