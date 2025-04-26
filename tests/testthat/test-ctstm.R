context("ctstm.R unit tests")
library("flexsurv")
library("mstate")
library("data.table")
rm(list = ls())

# Helper functions  ------------------------------------------------------------
pv <- function(z, r, t1, t2) {
  z * ((exp(-r * t1) - exp(-r * t2))/r)
}

# Strategies, population, and model structure  ---------------------------------
strategies <- data.table(strategy_id = c(1, 2))
patients <- data.table(patient_id = seq(1, 3),
                          age = c(45, 50, 60),
                          female = c(0, 0, 1))
states <- data.table(state_id = c(1, 2))
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients,
                        states = states)
n_samples <- 2

# Statistical models to estimate parameters  -----------------------------------
# Utility
utility_tbl <- stateval_tbl(data.frame(state_id = states$state_id,
                                       est = c(0.90, 0.55)),
                            dist = "fixed")

# Costs
## Medical
medcost_tbl <- stateval_tbl(data.frame(state_id = states$state_id,
                                       mean = c(800, 1500),
                                       se = c(100, 150)),
                            dist = "gamma")


## Drugs
drugcost_tbl <- stateval_tbl(tbl = data.frame(strategy_id = strategies$strategy_id,
                                           est = c(10000, 12500)),
                            dist = "fixed")

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
tmat_ebmt4 <- transMat(
  x = list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6), c(5, 6),
          c(), c()), 
  names = c("Tx", "Rec", "AE", "Rec+AE", "Rel", "Death")
)
msebmt <- msprep(data = ebmt4, trans = tmat_ebmt4, time = c(NA, "rec", "ae",
                 "recae", "rel", "srv"), status = c(NA, "rec.s", "ae.s", "recae.s",
                  "rel.s", "srv.s"), keep = c("match", "proph", "year", "agecl"))
covs <- c("match", "proph", "year", "agecl")
msebmt <- expand.covs(msebmt, covs, longnames = FALSE)
msfit <- flexsurvreg(Surv(time/365.25, status) ~ factor(trans), data = msebmt,
                           dist = "weibull")
msfit_cf <- flexsurvreg(Surv(Tstart/365.25, Tstop/365.25, status) ~ factor(trans),
                        data = msebmt,
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
  data.table(
    hesim:::C_ctstm_indiv_stateprobs(
      disprog, 
      t = time, 
      n_samples = 2,
      n_strategies = 2, 
      unique_strategy_id = unique(strategy_id),
      strategy_index = c(rep(0, 5), rep(1, 4), rep(0, 2), rep(1, 2)),
      n_grps = 1,
      unique_grp_id = 1,
      grp_index = rep(0, 13),
      n_states = 3,
      n_patients = 2
    )
  )
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

# Create transition model ------------------------------------------------------
# Transition specific models
msfit_list_data <- expand(hesim_dat)

test_that("create_IndivCtstmTrans has correct fields for transition specific models", {
  m <- create_IndivCtstmTrans(msfit_list, input_data = msfit_list_data, 
                              trans_mat = tmat_ebmt4, 
                              clock = "reset", n = n_samples)
  expect_equal(m$clock, "reset")
  expect_equal(m$trans_mat, tmat_ebmt4)
  expect_equal(m$params[[1]]$n_samples, n_samples)
  
  m <- create_IndivCtstmTrans(msfit_list, input_data = msfit_list_data, 
                              trans_mat = tmat_ebmt4, 
                              clock = "forward", n = n_samples)
  expect_equal(m$clock, "forward")
})

# Joint model
transitions <- create_trans_dt(tmat_ebmt4)
transitions[, trans := transition_id]
hesim_dat$transitions <- transitions
msfit_data <- expand(hesim_dat, by = c("strategies", "patients", "transitions"))

test_that("create_IndivCtstmTrans has correct fields for joint model", {
  m <- create_IndivCtstmTrans(msfit, input_data = msfit_data, 
                              trans_mat = tmat_ebmt4, 
                              clock = "reset", n = n_samples)
  expect_equal(m$clock, "reset")
  expect_equal(m$trans_mat, tmat_ebmt4)
  expect_equal(m$params$n_samples, n_samples)
  
  m <- create_IndivCtstmTrans(msfit_cf, input_data = msfit_data, 
                              trans_mat = tmat_ebmt4, 
                              clock = "forward", n = n_samples)
  expect_equal(m$clock, "forward")
})

# Simulate economic model ------------------------------------------------------
# Construct economic model
## Utility
utilmod <- create_StateVals(utility_tbl, n = n_samples, hesim_data = hesim_dat)

## Costs
medcostsmod <- create_StateVals(medcost_tbl, n = n_samples, hesim_data = hesim_dat)
drugcostsmod <- create_StateVals(drugcost_tbl, n = n_samples, hesim_data = hesim_dat) 

## Transitions
### With transition specific survival models
tmat <- rbind(c(NA, 1, 2),
              c(NA, NA, 3),
              c(NA, NA, NA))
test_that("IndivCtstmTrans - transition specific", {
  mstate_list <- create_IndivCtstmTrans(msfit_list, input_data = msfit_list_data, 
                                        trans_mat = tmat,
                                        uncertainty = "none")  
  
  # Hazard
  hesim_hazard <- mstate_list$hazard(3)
  expect_equal(hesim_hazard[transition_id == 1][1]$hazard,
               summary(msfit_list[[1]], type = "hazard", t = 3, ci = FALSE)[[1]][1, "est"])
  expect_equal(hesim_hazard[transition_id == 2][1]$hazard,
               summary(msfit_list[[2]], type = "hazard", t = 3, ci = FALSE)[[1]][1, "est"])

  # Cumulative hazard
  hesim_cumhazard <- mstate_list$cumhazard(5)
  expect_equal(hesim_cumhazard[transition_id == 1][1]$cumhazard,
               summary(msfit_list[[1]], type = "cumhaz", t = 5)[[1]][1, "est"])
  expect_equal(hesim_cumhazard[transition_id == 2][1]$cumhazard,
               summary(msfit_list[[2]], type = "cumhaz", t = 5)[[1]][1, "est"])
  
  # State probabilities only
  stprobs <- mstate_list$sim_stateprobs(t = c(0, 1, 2, 3))
  expect_true(inherits(stprobs, "data.table"))
  stprobs <- mstate_list$sim_stateprobs(t = c(0, 1, 2, 3),
                                        max_t = 2)
  expect_true(inherits(stprobs, "data.table"))
  
  # Disease progression and then state probabilities
  disprog <- mstate_list$sim_disease(max_t = 2)
  expect_true(inherits(disprog, "data.table"))
  stprobs <- mstate_list$sim_stateprobs(t = c(0, 2), disprog = disprog)
  expect_true(inherits(stprobs, "data.table"))
  
  # No errors
  expect_error(create_IndivCtstmTrans(msfit_list, input_data = msfit_list_data, 
                                      trans_mat = tmat,
                                      death_state = nrow(tmat),
                                      uncertainty = "none"),
               NA)
  
  # Errors
  expect_error(create_IndivCtstmTrans(msfit_list, input_data = msfit_list_data, 
                                      trans_mat = tmat,
                                      start_age = rep(55, nrow(dt_patients) + 1))) 
  expect_error(create_IndivCtstmTrans(msfit_list, input_data = msfit_list_data, 
                                      trans_mat = tmat,
                                      death_state = nrow(tmat) + 1,
                                      uncertainty = "none"))
  
  mstate_list2 <- mstate_list$clone()
  mstate_list2$trans_mat <- matrix(1)
  mstate_list2$trans_mat <-1
  expect_error(mstate_list2$sim_stateprobs(t = 2))
  mstate_list2$trans_mat <- matrix(seq(1, 6), nrow = 2)
  expect_error(mstate_list2$sim_stateprobs(t = 2))
})  

### With a joint model
test_that("IndivCtstmTrans - joint", {
  mstate <- create_IndivCtstmTrans(msfit, input_data = msfit_data, 
                                   trans_mat = tmat_ebmt4,
                                   uncertainty = "none")    
  # hazard
  hesim_hazard <- mstate$hazard(3)
  expect_equal(hesim_hazard[transition_id == 1][1]$hazard,
               summary(msfit, type = "hazard", t = 3, ci = FALSE)[[1]][1, "est"])
  expect_equal(hesim_hazard[transition_id == 4][1]$hazard,
               summary(msfit, type = "hazard", t = 3, ci = FALSE)[[4]][1, "est"])

  # Cumulative hazard
  hesim_cumhazard <- mstate$cumhazard(3)
  expect_equal(hesim_cumhazard[transition_id == 1][1]$cumhazard,
               summary(msfit, type = "cumhaz", t = 3, ci = FALSE)[[1]][1, "est"])
  expect_equal(hesim_cumhazard[transition_id == 4][1]$cumhazard,
               summary(msfit, type = "cumhaz", t = 3, ci = FALSE)[[4]][1, "est"])
  
  # Simulate state probabilities
  expect_true(inherits(mstate$sim_stateprobs(t = c(0, 1, 2)), "data.table")) 
})

# Simulate outcomes
## With transition specific survival models
test_that("Simulate disease progression with transition specific models", {
  mstate_list <- create_IndivCtstmTrans(msfit_list, input_data = msfit_list_data, 
                                         trans_mat = tmat,
                                         uncertainty = "none")
  ictstm <- IndivCtstm$new(trans_model = mstate_list,
                           utility_model = utilmod)
  
  # Errors
  expect_error(ictstm$sim_stateprobs())
  expect_error(ictstm$sim_disease(max_age = 12))
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
  mstate_list <- create_IndivCtstmTrans(msfit_list, input_data = msfit_list_data, 
                                         trans_mat = tmat,
                                         n = n_samples)  
  ictstm <- IndivCtstm$new(trans_model = mstate_list,
                           utility_model = utilmod,
                           cost_models = list(medical = medcostsmod, 
                                              drugs = drugcostsmod))
  ictstm$sim_disease()
  
  # Simulate QALYs
  ## Across patients
  ### Time constant utility
  #### Errors
  ##### Incorrect utility model
  ictstm2 <- ictstm$clone()
  ictstm2$utility_model <- 2
  expect_error(ictstm2$sim_qalys())
  
  ##### Cannot repeat discount rates
  expect_error(ictstm$sim_qalys(dr = c(.03, .03))$qalys)
  
  ##### Incorrect number of PSA samples
  ictstm2 <- ictstm$clone(deep = TRUE)
  ictstm2$utility_model$params$n_samples <- -1
  expect_error(
    ictstm2$sim_qalys(),
    paste0("The number of samples in each 'StateVals' model must equal the ",
           "number of samples in the 'disprog' object, which is 2.")
  )
  
  ##### Incorrect number of strategies
  ictstm2 <- ictstm$clone(deep = TRUE)
  ictstm2$utility_model$params$n_strategies <- -1
  expect_error(
    ictstm2$sim_qalys(),
    paste0("The number of strategies in each 'StateVals' model must equal the ",
           "number of strategies in the 'disprog' object, which is 2.")
  )
  
  ##### Incorrect number of patients
  ictstm2 <- ictstm$clone(deep = TRUE)
  ictstm2$utility_model$params$n_patients <- -1
  expect_error(
    ictstm2$sim_qalys(),
    paste0("The number of patients in each 'StateVals' model must equal the ",
           "number of patients in the 'disprog' object, which is 3.")
  )
  
  ##### Incorrect number of health states
  ictstm2 <- ictstm$clone(deep = TRUE)
  ictstm2$utility_model$params$n_states <- -1
  expect_error(
    ictstm2$sim_qalys(),
    paste0("The number of states in each 'StateVals' model must equal ",
           "the number of states in the 'disprog' object less the number ",
           "of absorbing states, which is 2."),
    fixed = TRUE
  )
  
  #### No errors
  expect_error(ictstm$sim_qalys(dr = .03)$qalys_, NA)
  
  ## By patient
  ### Time constant utility
  #### dr = .03
  ictstm$sim_qalys(dr = .03, by_patient = TRUE)$qalys_
  
  disprog1 <- ictstm$disprog_[sample == 1 & strategy_id == 1 & patient_id == 2]
  qalys1 <- ictstm$qalys_[sample == 1 & strategy_id == 1 & patient_id == 2]
  utilvals <- merge(data.frame(state_id = disprog1$from), utility_tbl, 
                    by = "state_id")
  qalys_expected <- sum(pv(utilvals$est, .03, disprog1$time_start, disprog1$time_stop))
  expect_equal(sum(qalys1$qalys), qalys_expected)
  
  #### dr = 0
  ictstm <- ictstm$clone(deep = TRUE)
  ictstm$utility_model$params$value <- matrix(1, 
                                              nrow = nrow(ictstm$utility_model$params$value),
                                              ncol = ncol(ictstm$utility_model$params$value))
  ictstm$sim_qalys(dr = 0, by_patient = TRUE)
  qalys <- ictstm$sim_qalys(dr = 0, by_patient = TRUE)$qalys_
  
  expect_equal(ictstm$disprog_[final == 1][sample == 1 & strategy_id == 2 & patient_id == 2, time_stop],
               sum(qalys[sample == 1 & strategy_id == 2 & patient_id == 2, qalys]))
  expect_true(all(qalys$qalys == qalys$lys)) # life-years are computed correctly
  
  ### Time varying utility
  t2 <- .2
  utility_tbl2 <- data.table(state_id = rep(states$state_id, 2),
                             time_start = c(0, 0, t2, t2),
                             est = c(.90, .55, .70, .35))
  utility_tbl2 <- stateval_tbl(utility_tbl2, dist = "fixed")
  utilmod2 <- create_StateVals(utility_tbl2, n = n_samples, hesim_data = hesim_dat)
  ictstm2 <- ictstm$clone()
  ictstm2$utility_model <- utilmod2
  ictstm2$sim_disease()
  ictstm2$sim_qalys(by_patient = TRUE, dr = .03)
  
  disprog1 <- ictstm2$disprog_[sample == 1 & strategy_id == 1 & patient_id == 2]
  qalys1 <- ictstm2$qalys_[sample == 1 & strategy_id == 1 & patient_id == 2]
  wlos1 <- rep(NA, nrow(disprog1))
  for (i in 1:2){
    if (i <= nrow(disprog1)){
      if (disprog1$time_start[i] < t2){
        z <- utility_tbl2[state_id == disprog1[i]$from & time_start == 0]$est
        wlos1[i] <- pv(z, .03, disprog1[i]$time_start, min(t2, disprog1[i]$time_stop))
        if (disprog1[i]$time_stop > t2){
          z <- utility_tbl2[state_id == disprog1[i]$from & time_start == t2]$est
          wlos1[i] <- wlos1[i] + pv(z, .03, t2, disprog1[i]$time_stop)
        }
      } else{
          z <- utility_tbl2[state_id == disprog1[i]$from & time_start == t2]$est
          wlos1[i] <- pv(z, .03, disprog1[i]$time_start, disprog1[i]$time_stop)
      }      
    } else{
      wlos1[i] <- 0
    }
  }
  expect_equal(wlos1, qalys1$qalys)

  # Simulate costs
  ## Errors
  ### Cost models must be lists of StateVal objects
  ictstm2$cost_models <- 2
  expect_error(ictstm2$sim_costs())
  ictstm2$cost_models <- list(2)
  expect_error(ictstm2$sim_costs())
  
  ## Working
  costs <- ictstm$sim_costs(dr = c(0, .03), max_t = c(Inf, 2))$costs_
  expect_true(inherits(costs, "data.table"))
  costs <- ictstm$sim_costs(dr = c(0, .03), max_t = c(Inf, Inf))$costs_
  expect_true(inherits(costs, "data.table"))
  costs <- ictstm$sim_costs(dr = c(0, .03), max_t = c(1, 0))$costs_
  expect_equal(costs[category == "drugs", costs], 
              rep(0, nrow(costs[category == "drugs"])),
              tolerance = .001, scale = 1)
  costs <- ictstm$sim_costs(dr = c(0, .03), max_t = c(0, 0))$costs_
  expect_equal(costs$costs, 
               rep(0, length(costs$costs)),
               tolerance = .001, scale = 1)
  
  costs <- ictstm$sim_costs(dr = c(0, .03), by_patient = TRUE)$costs_
  expect_equal(unique(costs$category), c("medical", "drugs"))
  expect_equal(unique(costs$dr), c(0, .03))
  
  ## Time-varying costs with reset
  ### Simulate
  drugcost_tbl_tv <- stateval_tbl(tbl = data.frame(strategy_id = rep(strategies$strategy_id, each = 2),
                                                 time_start = c(0, 1.2, 0, 1.2),
                                                  est = c(10000, 500, 12500, 250)),
                                  dist = "fixed")
  drugcostsmod_tv <- create_StateVals(drugcost_tbl_tv, n = n_samples, time_reset = TRUE,
                                      hesim_data = hesim_dat) 
  ictstm2$cost_models <- list(medical = medcostsmod, 
                              drugs = drugcostsmod_tv)
  ictstm2$sim_costs(dr = 0, by_patient = TRUE)
  ictstm2$disprog_[, time_elapsed := time_stop - time_start]
  
  ### Test
  test_tv_cost <- function(state){
    row <- ictstm2$disprog_[from == state][1]
    drug_costs <- drugcost_tbl_tv[strategy_id == row$strategy_id]
    expected_costs <- 0
    j <- 1
    while(j <= nrow(drug_costs) & (drug_costs[j]$time_start < row$time_elapsed)){
      time_j <- min(row$time_elapsed - drug_costs[j]$time_start, 
                    drug_costs[j]$time_stop - drug_costs[j]$time_start)
      expected_costs <- expected_costs + time_j * drug_costs[j, est]
      j <- j + 1
    }
    expect_equal(ictstm2$costs_[category == "drugs" & sample == row$sample & 
                                strategy_id == row$strategy_id & 
                                patient_id == row$patient_id & 
                                state_id == row$from]$costs,
               expected_costs)     
  }
  test_tv_cost(1)
  test_tv_cost(2)
  
  ## Using method = "starting"
  drugcost_tbl_tv2 <- stateval_tbl(tbl = data.frame(strategy_id = rep(strategies$strategy_id, 
                                                                      each = 3),
                                                   time_start = c(0, 1.2, 4, 
                                                                  0, 1.2, 4),
                                                   est = c(10000, 500, 0, 
                                                           12500, 250, 0)),
                                  dist = "fixed")
  drugcostsmod_tv2 <- create_StateVals(drugcost_tbl_tv2, n = n_samples, 
                                       time_reset = FALSE, method = "starting",
                                       hesim_data = hesim_dat) 
  ictstm2$cost_models$drugs <- drugcostsmod_tv2
  ictstm2$sim_costs(dr = .03, by_patient = TRUE)
  
  test_starting <- function(model){
    cost_params <- model$cost_models$drugs$params 
    n_obs <- nrow(cost_params$value)
    n_samples <- ncol(cost_params$value)
    cost_tbl <- data.table(value = c(cost_params$value),
                           sample = rep(1:n_samples, 
                                        each = n_obs),
                           strategy_id = rep(cost_params$strategy_id,
                                             times = n_samples),
                           patient_id = rep(cost_params$patient_id,
                                            times = n_samples),
                           time_id = rep(cost_params$time_id,
                                         times = n_samples),
                           from = rep(cost_params$state_id,
                                      times = n_samples))
    new_disprog <- copy(model$disprog_)
    new_disprog[, time_id := findInterval(time_start, 
                                          cost_params$time_intervals$time_start)]
    new_disprog <- merge(new_disprog, cost_tbl,
                         by = c("sample", "strategy_id", "patient_id", "from",
                                "time_id"))
    new_disprog[, R_cost := exp(-.03 * time_start) * value]
    setnames(new_disprog, "from", "state_id")
    R_costs <- new_disprog[, .(R_costs = sum(R_cost)), 
                           by  = c("sample", "strategy_id", "patient_id", "state_id")]
    model_costs <- model$costs_[category == "drugs"]
    costs <- merge(model_costs, R_costs, 
                   by = c("sample", "strategy_id", "patient_id", "state_id"),
                   all.x = TRUE)
    costs[is.na(R_costs), R_costs := 0]
    expect_equal(costs$costs, costs$R_costs)
  }
  test_starting(ictstm2)
   

  ## Using method = "transition"
  drugcost_tbl2 <- stateval_tbl(tbl = data.frame(state_id = 1:3,
                                                 est = c(10, 20, 30)),
                                dist = "fixed")
  drugcostsmod2 <- create_StateVals(drugcost_tbl2, n = n_samples,
                                    time_reset = FALSE, method = "transition",
                                    hesim_data = hesim_dat) 
  ictstm2 <- ictstm$clone()
  ictstm2$cost_models <- list(drugs=drugcostsmod2)
  ictstm2$sim_costs(dr = 0, by_patient = TRUE)
  test_transition_costs <- function(model){
      model$disprog_[, costs := ifelse(from==1 & to==2, 10,
                                ifelse(from==1 & to==3, 20, 30))]
      expect_equal(model$disprog_[, sum(costs), by=from]$V1,
                   model$costs_[, sum(costs), by=state_id]$V1)
  }
  test_transition_costs(ictstm2)
  
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
  
  # Cost-effectiveness analysis
  cea <- cea(ce_summary, dr_costs = 0, dr_qalys = 0)
  expect_true("mce" %in% names(cea))
  cea <- cea(ce_summary, dr_qalys = 0, dr_costs = .03)
  expect_true("mce" %in% names(cea))
  
  cea_pw <- cea_pw(ce_summary, comparator = 1, 
                   dr_costs = 0, dr_qalys = 0)
  expect_true("ceac" %in% names(cea_pw))
})

## With a joint survival model
test_that("IndivCtstm - joint", {
  # Clock-reset
  mstate <- create_IndivCtstmTrans(msfit, input_data = msfit_data, 
                                   trans_mat = tmat_ebmt4,
                                   n = n_samples)      
  ictstm <- IndivCtstm$new(trans_model = mstate)
  
  # Simulate disease progression
  expect_error(ictstm$sim_disease()$disprog_, NA)
  disprog <- ictstm$sim_disease()$disprog_
  expect_error(ictstm$sim_disease(progress = 1), NA) # outputs sample = 1 and sample = 2
  expect_true(is.data.table(disprog))
  
  # Simulate state probabilities
  stprobs <- ictstm$sim_stateprobs(t = c(0, 1, 2, 3))$stateprobs_
  expect_true(all(stprobs[t == 0 & state_id == 1, prob] ==  1))
  expect_true(all(stprobs[t == 0 & state_id != 1, prob] ==  0))
  
  # Clock-forward
  mstate_cf <- create_IndivCtstmTrans(msfit_cf, input_data = msfit_data, 
                                      trans_mat = tmat_ebmt4,
                                      clock = "forward", n = n_samples)
  ictstm <- IndivCtstm$new(trans_model = mstate_cf)
  disprog <- ictstm$sim_disease()$disprog_
  expect_true(is.data.table(disprog))
  
  # Mixture - states
  params <- create_params(msfit, n = 2)
  msfit_data <- cbind(msfit_data, model.matrix(~factor(trans), msfit_data))
  mstate_mix <- create_IndivCtstmTrans(params, input_data = msfit_data, 
                                       trans_mat = tmat_ebmt4,
                                       clock = "mix", reset_states = 1:nrow(tmat_ebmt4))  
  expect_true(inherits(mstate_mix, "IndivCtstmTrans"))
  ictstm <- IndivCtstm$new(trans_model = mstate_mix)
  disprog <- ictstm$sim_disease()$disprog_
  expect_true(is.data.table(disprog))

  # Mixture - transitions
  test_mixt <- function(trans_type) {
    mstate_mixt <- create_IndivCtstmTrans(params, input_data = msfit_data,
                                          trans_mat = tmat_ebmt4,
                                          clock = "mixt",
                                          transition_types = trans_type)
    expect_true(inherits(mstate_mixt, "IndivCtstmTrans"))
    ictstm <- IndivCtstm$new(trans_model = mstate_mixt)
    disprog <- ictstm$sim_disease()$disprog_
    expect_true(is.data.table(disprog))
  }
  test_mixt(rep("reset", max(tmat_ebmt4, na.rm=TRUE)))
  test_mixt(rep("age", max(tmat_ebmt4, na.rm=TRUE)))
  test_mixt(rep("time", max(tmat_ebmt4, na.rm=TRUE)))
})

## With fractional polynomial or survival spline from parameters object
### Input data
edat <- copy(msfit_data[patient_id == 1])
edat$intercept <- 1

### Simulate disease progression
sim_disprog <- function(params){
  mstate_fp <- create_IndivCtstmTrans(params, input_data = edat, 
                                      trans_mat = tmat_ebmt4)
  ictstm <- IndivCtstm$new(trans_model = mstate_fp)
  return(ictstm$sim_disease()$disprog_)
}
  
test_that("IndivCtstm - fracpoly", {
  # Initial parameters
  coefs_fp <- list(gamma0 = matrix(log(1/5), nrow = 2, ncol = 1),
                   gamma1 = matrix(1, nrow = 2, ncol = 1))
  colnames(coefs_fp$gamma0) <- colnames(coefs_fp$gamma1) <- c("intercept")
  fp_args <- list(coefs = coefs_fp, dist = "fracpoly", aux = list(powers = 0))
  params_fp <- do.call("params_surv", fp_args)
  
  # Working
  ## Inverse CDF and quadrature
  expect_true(is.data.table(sim_disprog(params_fp)))
  
  ## Inverse CDF and riemann
  params_fp$aux$cumhaz_method <- "riemann"
  params_fp$aux$step <- 1/12
  expect_true(is.data.table(sim_disprog(params_fp)))
  
  ## Sample and riemann
  params_fp$aux$random_method <- "discrete"
  expect_true(is.data.table(sim_disprog(params_fp)))
  
  ## Sample and quadrature
  params_fp$aux$cumhaz_method <- "quad"
  expect_true(is.data.table(sim_disprog(params_fp)))
  
  # Errors
  ## Need step size
  ### (1)
  fp_args$aux$random_method <- "discrete"
  expect_error(do.call("params_surv", fp_args))
  
  ### (2)
  fp_args$aux$random_method <- "invcdf"
  fp_args$aux$cumhaz_method <- "riemann"
  expect_error(do.call("params_surv", fp_args))
})

## With survival splines from parameters object
test_that("IndivCtstm - survspline", {
  # Initial parameters
  coefs_spline <- list(gamma0 = matrix(-2, nrow = 2, ncol = 1),
                       gamma1 = matrix(1, nrow = 2, ncol = 1))  
  colnames(coefs_spline$gamma0) <- colnames(coefs_spline$gamma1) <- c("intercept")
  spline_args <- list(coefs = coefs_spline,
                      dist = "survspline",
                      aux = list(knots = c(-10, 10),
                                 scale = "log_hazard",
                                 timescale = "log",
                                 random_method = "invcdf"))
  params_spline <- do.call("params_surv", spline_args)
  
  # Working
  expect_equal(params_spline$aux$random_method, "invcdf")
  expect_equal(params_spline$aux$cumhaz_method, "quad")
  
  # Warnings
  spline_args$aux$random_method <- "sample"
  spline_args$aux$scale <- "log_cumhazard"
  expect_warning(
    do.call("params_surv", spline_args),
    "'random_method' = 'sample' is deprecated. Use 'discrete' instead."
  )
  
  # Errors
  ## Need step size
  ### (1)
  spline_args$aux$random_method <- "discrete"
  spline_args$aux$scale <- "log_hazard"
  expect_error(
    do.call("params_surv", spline_args),
    "'step' must be specified"
  )
  
  ### (2)
  spline_args$aux$random_method <- "invcdf"
  spline_args$aux$cumhaz_method <- "riemann"
  expect_error(
    do.call("params_surv", spline_args),
    "'step' must be specified"
  ) 
})
