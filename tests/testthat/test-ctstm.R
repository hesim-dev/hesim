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
utility_dist <- matrix(rbeta(n_samples * 2, shape1 = beta_params$shape1,
                                  shape2 = beta_params$shape2),
                            nrow = n_samples, byrow = TRUE)

# Costs
## Medical
gamma_params <- mom_gamma(ctstm3_exdata$costs$medical$mean, 
                        ctstm3_exdata$costs$medical$se)
medcosts_dist <- matrix(rgamma(n_samples * 2, shape = gamma_params$shape,
                                scale = gamma_params$scale),
                            nrow = n_samples, byrow = TRUE)

## Drugs
drugcosts <- ctstm3_exdata$costs$drugs$costs

# Clock-reset multi-state model
## Separate 
dat <- data.table(bosms3)

### fits
msfit_list <- vector(length = 3, mode = "list")
for (i in 1:length(msfit_list)){
  msfit_list[[i]] <- flexsurvreg(Surv(years, status) ~ 1, data = dat[trans == i],
                           dist = "exp")
}

msfit_list_wei <- vector(length = 3, mode = "list")
for (i in 1:length(msfit_list)){
  msfit_list_wei[[i]] <- flexsurvreg(Surv(years, status) ~ 1, data = dat[trans == i],
                           dist = "weibull")
}

### Flexsurvreg lists
msfit_list_mix <- flexsurvreg_list(msfit_list[[1]], msfit_list[[2]],
                                   msfit_list_wei[[3]])
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
sample <-  c(rep(1, 9), rep(2, 4)) 
strategy_id <- c(rep(2, 5), rep(3, 4), rep(2, 2), rep(3, 2))
patient_id <- c(rep(2, 3), rep(3, 2), rep(2, 2), rep(3, 2), rep(2, 2), rep(3, 2)) 
from <- c(1, 2, 1, 1, 2, 1, 2, 2, 1, rep(1, 4)) 
to <- c(2, 1, 3, 2, 3, 2, 3, 1, 3, rep(3, 4)) 
final <- c(0, 0, 1, 0, 1, 0, 1, 0, 1, rep(1, 4))
time_start <- c(0, 2, 5, 0, 3, 2, 4, 0, 5, rep(0, 4))
time_stop <- c(2, 5, 7, 3, 8, 4, 8, 5, 11, rep(2, 2), rep(2.6, 2))

disprog <- data.frame(sample = sample, 
                      strategy_id = strategy_id,
                      patient_id = patient_id,
                      from = from,
                      to = to,
                      final = final,
                      time_start = time_start, 
                      time_stop = time_stop)

stprob_fun <- function(time){
  data.table(hesim:::C_ctstm_indiv_stateprobs(disprog, t = time, n_samples = 2,
                                              n_strategies = 2, 
                                              unique_strategy_id = unique(strategy_id),
                                              strategy_index = c(rep(0, 5), rep(1, 4), rep(0, 2), rep(1, 2)),
                                              n_states = 3,
                                              n_patients = 2))
}


test_that("C_ctstm_indiv_stateprobs", {
  # Time 2.5
  stprobs <- stprob_fun(2.5)
  expect_equal(stprobs[sample == 0 & state_id == 1 & strategy_id == 2, prob], .5)
  expect_equal(stprobs[sample == 0 & state_id == 2 & strategy_id == 2, prob], 0)
  expect_equal(stprobs[sample == 0 & state_id == 0 & strategy_id == 3, prob], .5)
  expect_equal(stprobs[sample == 0 & state_id == 1 & strategy_id == 3, prob], .5)
  expect_equal(stprobs[sample == 1 & state_id == 2 & strategy_id == 2, prob], 1)
  expect_equal(stprobs[sample == 1 & state_id == 0 & strategy_id == 3, prob], 1)
  
  # Time 3
  stprobs <- stprob_fun(3)
  expect_equal(stprobs[sample == 0 & state_id == 1 & strategy_id == 2, prob], 1)
  expect_equal(stprobs[sample == 0 & state_id == 1 & strategy_id == 3, prob], .5)
  
  # Time 10
  stprobs <- stprob_fun(10)
  expect_equal(stprobs[sample == 0 & state_id == 1 & strategy_id == 2, prob], 0)
  expect_equal(stprobs[sample == 0 & state_id == 2 & strategy_id == 2, prob], 1)
  expect_equal(stprobs[sample == 0 & state_id == 0 & strategy_id == 3, prob], .5)
  expect_equal(stprobs[sample == 0 & state_id == 2 & strategy_id == 3, prob], .5)
})

# Simulate economic model ------------------------------------------------------
# Construct a model structure 
## Simulation data
dt_strategies <- data.table(strategy_id = c(1, 2))
dt_patients <- data.table(patient_id = seq(1, 3),
                          age = c(45, 50, 60),
                          female = c(0, 0, 1))
dt_states <- data.table(state_id = c(1, 2))
hesim_dat <- hesim_data(strategies = dt_strategies,
                        patients = dt_patients,
                        states = dt_states)

## Utility and cost models
### Utility
utility_means <- stateval_means(values = utility_dist,
                                strategy_id = dt_strategies$strategy_id,
                                patient_id = dt_patients$patient_id)
utilmod <- create_StateVals(utility_means)

### Costs
#### Medical
medcosts_means <- stateval_means(values = medcosts_dist,
                                strategy_id = dt_strategies$strategy_id,
                                patient_id = dt_patients$patient_id)
medcostsmod <- create_StateVals(medcosts_means)

#### Drug
drugcost_array <- array(NA, dim = c(n_samples, 2, 2))
drugcost_array[, , 1] <- drugcosts[1]
drugcost_array[, , 2] <- drugcosts[2]
drugcost_means <- stateval_means(values = drugcost_array,
                                strategy_id = dt_strategies$strategy_id,
                                patient_id = dt_patients$patient_id)
drugcostsmod <- create_StateVals(drugcost_means)

## The health state transitions
### With transition specific survival models
msfit_list_data <- expand(hesim_dat)
tmat <- rbind(c(NA, 1, 2),
              c(NA, NA, 3),
              c(NA, NA, NA))

test_that("IndivCtstmTrans - transition specific", {
  mstate_list <- create_IndivCtstmTrans(msfit_list, data = msfit_list_data, trans_mat = tmat,
                                      point_estimate = TRUE)  
  
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
  stprobs <- mstate_list$sim_stateprobs(t = c(0, 1, 2, 3))
  expect_true(inherits(stprobs, "data.table"))
  stprobs <- mstate_list$sim_stateprobs(t = c(0, 1, 2, 3),
                                        max_t = 2)
  expect_true(inherits(stprobs, "data.table"))
  
  # No errors
  expect_error(create_IndivCtstmTrans(msfit_list, data = msfit_list_data, trans_mat = tmat,
                                death_state = nrow(tmat),
                                point_estimate = TRUE),
               NA)
  
  # Errors
  expect_error(create_IndivCtstmTrans(msfit_list, data = msfit_list_data, 
                                      trans_mat = tmat,
                                      start_age = rep(55, nrow(dt_patients) + 1))) 
  expect_error(create_IndivCtstmTrans(msfit_list, data = msfit_list_data, trans_mat = tmat,
                                death_state = nrow(tmat) + 1,
                                point_estimate = TRUE))
  
  mstate_list2 <- mstate_list$clone()
  mstate_list2$trans_mat <- matrix(1)
  expect_error(mstate_list2$sim_disease(t = 2))
  mstate_list2$trans_mat <-1
  expect_error(mstate_list2$sim_stateprobs(t = 2))
  mstate_list2$trans_mat <- matrix(seq(1, 6), nrow = 2)
  expect_error(mstate_list2$sim_stateprobs(t = 2))
})  

### With a joint model
dt_transitions <- create_trans_dt(tmat_ebmt4)
dt_transitions[, trans := transition_id]
hesim_dat$transitions <- dt_transitions
msfit_data <- expand(hesim_dat, by = c("strategies", "patients", "transitions"))

test_that("create_IndivCtstmTrans - joint", {
  # From parameter object
  params <- create_params(msfit, n = 2)
  tmp_data <- msfit_data
  msfit_data <- cbind(msfit_data, model.matrix(~factor(trans), msfit_data), shape = 1, scale = 1)
  obj <- create_IndivCtstmTrans(params, data = msfit_data, trans_mat = tmat_ebmt4)
  expect_true(inherits(obj, "IndivCtstmTrans"))
})

test_that("IndivCtstmTrans - joint", {
  mstate <- create_IndivCtstmTrans(msfit, data = msfit_data, trans_mat = tmat_ebmt4,
                            point_estimate = TRUE)    
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
  
  # Simulate state probabilities
  expect_true(inherits(mstate$sim_stateprobs(t = c(0, 1, 2)), "data.table")) 
})

# Simulate outcomes
## With transition specific survival models
test_that("Simulate disease progression with transition specific models", {
  mstate_list <- create_IndivCtstmTrans(msfit_list, data = msfit_list_data, 
                                         trans_mat = tmat,
                                         point_estimate = TRUE)
  ictstm <- IndivCtstm$new(trans_model = mstate_list,
                           utility_model = utilmod)
  
  # Errors
  expect_error(ictstm$sim_stateprobs())
  ictstm2 <- IndivCtstm$new(trans_model = 2)
  expect_error(ictstm2$sim_disease())
  
  # Base case simulation
  ## Simulate disease progression
  set.seed(101)
  disprog <- ictstm$sim_disease()$disprog_
  expect_true(inherits(disprog, "data.table"))
  
  ### Time from first state
  set.seed(101)
  time_1a <- rexp(1, rate = msfit_list[[1]]$res[, "est"])
  time_1b <- rexp(1, rate = msfit_list[[2]]$res[, "est"])
  time1 <- min(time_1a, time_1b)
  state1 <- which.min(c(time_1a, time_1b)) + 1
  expect_equal(disprog[1, time_stop], time1)
  expect_equal(disprog[1, to], state1)
  
  #### Time from second state
  time2 <- rexp(1, rate = msfit_list[[3]]$res[, "est"])
  expect_equal(disprog[2, time_stop],  time1 + time2)  
  
  ## Simulate state probabilities
  stprobs <- ictstm$sim_stateprobs(t = c(0, 1, 2, 3))$stateprobs_
  expect_true(all(stprobs[t == 0 & state_id == 1, prob] ==  1))
  expect_true(all(stprobs[t == 0 & state_id != 1, prob] ==  0))
  
  # Simulation with other cases 
  ## Maximum time = 2
  disprog <- ictstm$sim_disease(max_t = 2)$disprog_
  expect_equal(ictstm$stateprobs_, NULL)
  expect_equal(max(disprog$time_stop), 2)
  expect_true(all(disprog[final == 1 & time_stop == 2 & from == 1, to] == 1)) # All should have remained in initial state
  expect_true(all(disprog[final == 1 & time_stop == 2 & from == 2, to] == 2))
  
  ## Maximum age = 43 (i.e., max_t = 5)
  disprog <- ictstm$sim_disease(max_age = 43)$disprog_
  expect_true(all(disprog[time_stop == 5 & final == 1, to] == 3)) # All should have moved to death state
})

test_that("Simulate costs and QALYs", {
  mstate_list <- create_IndivCtstmTrans(msfit_list, data = msfit_list_data, 
                                         trans_mat = tmat,
                                         n = n_samples)  
  ictstm <- IndivCtstm$new(trans_model = mstate_list,
                           utility_model = utilmod,
                           cost_models = list(medical = medcostsmod, 
                                              drugs = drugcostsmod))
  ictstm$sim_disease()
  
  # Simulate QALYs
  ## Errors
  ictstm2 <- ictstm$clone()
  ictstm2$utility_model <- 2
  expect_error(ictstm2$sim_qalys())
  
  ## No errors
  expect_error(ictstm$sim_qalys(dr = .03)$qalys_, NA)
  
  ## By patient
  ### dr = .03
  ictstm$sim_qalys(dr = .03, by_patient = TRUE)$qalys_
  
  disprog1 <- ictstm$disprog_[sample == 1 & strategy_id == 1 & patient_id == 2]
  qalys1 <- ictstm$qalys_[sample == 1 & strategy_id == 1 & patient_id == 2]
  utilvals <- utility_means$values[1, disprog1$from] 
  qalys_expected <- sum(pv(utilvals, .03, disprog1$time_start, disprog1$time_stop))
  expect_equal(sum(qalys1$qalys), qalys_expected)
  
  ### dr = 0
  ictstm <- ictstm$clone(deep = TRUE)
  ictstm$utility_model$params$mu <- matrix(1, nrow = nrow(ictstm$utility_model$params$mu),
                                           ncol = ncol(ictstm$utility_model$params$mu))
  ictstm$sim_qalys(dr = 0, by_patient = TRUE)
  qalys <- ictstm$sim_qalys(dr = 0, by_patient = TRUE)$qalys_
  
  expect_equal(ictstm$disprog_[final == 1][sample == 1 & strategy_id == 2 & patient_id == 2, time_stop],
               sum(qalys[sample == 1 & strategy_id == 2 & patient_id == 2, qalys]))
  expect_true(all(qalys$qalys == qalys$lys)) # life-years are computed correctly
  
  # Simulate costs
  # Errors
  ictstm2$cost_models <- 2
  expect_error(ictstm2$sim_costs())
  ictstm2$cost_models <- list(2)
  expect_error(ictstm2$sim_costs())
  
  # Working
  costs <- ictstm$sim_costs(dr = c(0, .03), max_t = c(Inf, 2))$costs_
  expect_true(inherits(costs, "data.table"))
  costs <- ictstm$sim_costs(dr = c(0, .03), max_t = c(Inf, Inf))$costs_
  expect_true(inherits(costs, "data.table"))
  costs <- ictstm$sim_costs(dr = c(0, .03), max_t = c(1, 0))$costs_
  expect_true(all(costs[category == "drugs", costs] == 0))
  costs <- ictstm$sim_costs(dr = c(0, .03), max_t = c(0, 0))$costs_
  expect_true(all(costs$costs == 0))
  
  costs <- ictstm$sim_costs(dr = c(0, .03), by_patient = TRUE)$costs_
  expect_equal(unique(costs$category), c("medical", "drugs"))
  expect_equal(unique(costs$dr), c(0, .03))
   
  # Summarize costs and QALYs
  ## By patient = TRUE
  ce_summary <- ictstm$summarize()
  expect_true(inherits(ce_summary, "ce"))
  ictstm$sim_disease()
  expect_error(ictstm$summarize())
  ictstm$sim_costs()
  expect_error(ictstm$summarize())
  
  ## By patient = FALSE
  ictstm$sim_qalys(by_patient = TRUE)
  ictstm$sim_costs(by_patient = TRUE)
  expect_error(ictstm$summarize(), NA)
})

## With a joint survival model
test_that("IndivCtstm", {
  mstate <- create_IndivCtstmTrans(msfit, data = msfit_data, trans_mat = tmat_ebmt4,
                                    n = n_samples)      
  ictstm <- IndivCtstm$new(trans_model = mstate)
  
  # Simulate disease progression
  expect_error(ictstm$sim_disease()$disprog_, NA)
  disprog <- ictstm$sim_disease()$disprog_
  expect_true(is.data.table(disprog))
  
  # Simulate state probabilities
  stprobs <- ictstm$sim_stateprobs(t = c(0, 1, 2, 3))$stateprobs_
  expect_true(all(stprobs[t == 0 & state_id == 1, prob] ==  1))
  expect_true(all(stprobs[t == 0 & state_id != 1, prob] ==  0))
})


