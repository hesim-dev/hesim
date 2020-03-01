context("statevals.R unit tests")
library("data.table")
rm(list = ls())

# Simulation data
strategies <- data.table(strategy_id = c(1, 2))
n_strategies <- nrow(strategies)
patients <- data.table(patient_id = seq(1, 3),
                       grp = c(1, 1, 2),
                       age = c(45, 50, 60),
                       female = c(0, 0, 1))
n_patients <- nrow(patients)
states <- data.frame(state_id =  seq(1, 3),
                     state_name = paste0("state", seq(1, 3)))
n_states <- nrow(states)
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients,
                        states = states)
N <- 5

# stateval_tbl -----------------------------------------------------------------
n_grps <- 2
strategy_id <- rep(strategies$strategy_id, each = n_grps * n_states)
state_id <- rep(states$state_id, times = n_grps * n_strategies)
grp <- rep(rep(1:n_grps, each = n_states), times = n_strategies)
mean <- c(.80, .70, .60,    # Strategy 1, Group 1
          .70, .60, .50,    # Strategy 1, Group 2
          .90, .80, .70,    # Strategy 2, Group 1
          .85, .75, .65)    # Strategy 2, Group 2
se <- c(.20, .15, .10,
        .15, .10, .05,
        .23, .18, .13,
        .25, .20, .15)
beta_params <- mom_beta(mean, se)
shape1 <- beta_params$shape1
shape2 <- beta_params$shape2
tbl <- data.table(strategy_id = strategy_id, grp = grp,
                  state_id = state_id, mean = mean, se = se,
                  sd = se,
                  shape1 = shape1,
                  shape2 = shape2, min = .3, max = .5)

test_that("stateval_tbl with state_id and built-in distributions", {
  # Beta distribution
  stateval_tbl <- stateval_tbl(tbl[strategy_id == 1 & grp == 1,
                                   .(state_id, mean, se)], 
                               dist = "beta", 
                               hesim_data = hesim_dat) 
  mod <- create_StateVals(stateval_tbl, n = 2)
  expect_equal(nrow(mod$params$value), nrow(stateval_tbl) * 3 * 2)
  expect_equal(ncol(mod$params$value), 2)
  expect_true(all(mod$params$value[c(1, 4, 7, 10), 1] == mod$params$value[1, 1]))
  
  ## time_reset = TRUE
  mod <- create_StateVals(stateval_tbl, n = 2, time_reset = TRUE)
  expect_true(mod$time_reset)
  
  # Uniform distribution
  stateval_tbl <- stateval_tbl(tbl[strategy_id == 1 & grp == 1,
                                   .(state_id, min, max)], 
                               dist = "unif",  
                               hesim_data = hesim_dat)
  mod <- create_StateVals(stateval_tbl, n = 2)
  expect_true(all(mod$params$value >= .3 & mod$params$value <= .5))
  
  # Normal distribution
  stateval_tbl <- stateval_tbl(tbl[strategy_id == 1 & grp == 1,
                                   .(state_id, mean, sd)], 
                               dist = "norm",  
                               hesim_data = hesim_dat)  
  mod <- create_StateVals(stateval_tbl, n = 2) 
  expect_equal(nrow(mod$params$value), 
               nrow(stateval_tbl) * n_patients * n_strategies) 
})
  
test_that("stateval_tbl with state_id and custom distribution", {
  tbl2 <- data.table(state_id = rep(states$state_id, 2),
                     sample = rep(c(1, 2), each = n_states),
                     value = rnorm(6, 4))
  expect_error(stateval_tbl(tbl2, 
                            hesim_data = hesim_dat),
               "If 'sample' is in 'tbl', then 'dist' must equal 'custom'")
  stateval_tbl <- stateval_tbl(tbl2, 
                               dist = "custom",
                               hesim_data = hesim_dat)
  
  # n is less than the number of samples in tbl
  expect_warning(mod <- create_StateVals(stateval_tbl, n = 3))
  expect_equal(ncol(mod$params$value), 3)
  expect_true(mod$params$value[1, 1] == tbl2[state_id == 1 & sample == 1, value] |
                mod$params$value[1, 1] == tbl2[state_id == 1 & sample == 2, value]) 
  
  # n equals number of samples in tbl
  mod <- create_StateVals(stateval_tbl, n = 2) 
  expect_equal(ncol(mod$params$value), 2)
  expect_equal(nrow(mod$params$value), 
               n_states * n_patients * n_strategies) 
  expect_equal(mod$params$value[5, 1],
               tbl2[state_id == 2 & sample == 1, value])
  
  # n is less than the number of samples in tbl
  mod <- create_StateVals(stateval_tbl, n = 1) 
  expect_equal(ncol(mod$params$value), 1)
  expect_true(mod$params$value[3, 1] == tbl2[state_id == 3 & sample == 1, value] |
                mod$params$value[3, 1] == tbl2[state_id == 3 & sample == 2, value]) 
})

test_that("stateval_tbl with strategy_id", {
  # Gamma distribution
  stateval_tbl <- stateval_tbl(tbl[state_id == 3 & grp == 1,
                                   .(strategy_id, mean, se)], 
                               dist = "gamma",
                               hesim_data = hesim_dat)   
  mod <- create_StateVals(stateval_tbl, n = 2)
  expect_equal(nrow(mod$params$value), nrow(stateval_tbl) * 3 * 3)
  expect_equal(mod$params$state_id, rep(c(1, 2, 3), 3 * 2))
  expect_true(all(mod$params$value[1:9, 2] == mod$params$value[1, 2]))
})

test_that("stateval_tbl with patient_id", {
  tbl[, patient_id := grp]
  stateval_tbl <- stateval_tbl(tbl[state_id == 3 & grp == 1,
                                   .(strategy_id, mean, se)], 
                               dist = "gamma",
                               hesim_data = hesim_dat)   
  mod <- create_StateVals(stateval_tbl, n = 2)
  expect_equal(nrow(mod$params$value), nrow(stateval_tbl) * 3 * 3)
  expect_equal(mod$params$state_id, rep(c(1, 2, 3), 3 * 2))
  expect_true(all(mod$params$value[1:9, 2] == mod$params$value[1, 2]))
})

test_that("stateval_tbl with grp_var", {
  stateval_tbl <- stateval_tbl(tbl[, .(strategy_id, grp, state_id,
                                       shape1, shape2)], 
                               dist = "beta", 
                               hesim_data = hesim_dat,
                               grp_var = "grp")  
  mod <- create_StateVals(stateval_tbl, n = 1)
  expect_equal(ncol(mod$params$value), 1)
  expect_equal(mod$params$patient_id,
               rep(rep(patients$patient_id, each = nrow(states)), 
                   nrow(strategies))) 
})
  
test_that("stateval_tbl with state_id and time", {
  tbl2 <- tbl[strategy_id == 1]
  tbl2[, time_start := ifelse(grp == 1, 0, 4)]
  tbl2[, time_stop := ifelse(grp == 1, 4, 10)]
  stateval_tbl <- stateval_tbl(tbl2[,.(state_id, time_start, time_stop,
                                       mean, se)], 
                               dist = "beta",
                               hesim_data = hesim_dat)
  mod <- create_StateVals(stateval_tbl, n = 2)
  expect_equal(nrow(mod$params$value),
               n_states * n_strategies * n_patients * 2)
})

test_that("Errors wuth stateval_tbl", {
  # Don't have required tables in hesim_data
  expect_error(stateval_tbl(tbl[,.(strategy_id, grp, mean)], 
                            dist = "beta"),
               paste0("If 'state_id' is not a column in 'tbl', then ",
                       "'hesim_data' must be included as an argument and ", 
                       "'states' must be an element of 'hesim_data'."))
  expect_error(stateval_tbl(tbl[,.(state_id, grp, mean)], 
                            dist = "beta"),
               paste0("If 'strategy_id' is not a column in 'tbl', then ",
                      "'hesim_data' must be included as an argument and ", 
                      "'strategies' must be an element of 'hesim_data'."))  
  expect_error(stateval_tbl(tbl[,.(state_id, strategy_id, grp, mean)], 
                            dist = "beta"),
               paste0("If 'patient_id' is not a column in 'tbl', then ",
                      "'hesim_data' must be included as an argument and ", 
                      "'patients' must be an element of 'hesim_data'."))  
  
  # Correct number of rows
  row_msg <- paste0("There must only be one row for each combination ",
                    "of strategy_id and state_id in 'tbl'.")
  tbl2 <- tbl[,.(strategy_id, state_id, grp, mean, se)]
  tbl2 <- rbind(tbl2, data.table(strategy_id = 2, state_id = 1, grp = 1, 
                                 mean = .5, se = .1))
  expect_error(stateval_tbl(tbl2, dist = "beta",
                            hesim_data = hesim_dat),
               row_msg)
  
  tbl2 <- tbl2[(strategy_id == 1 & state_id == 1) |
                 (strategy_id == 2 & state_id == 2)]
  expect_error(stateval_tbl(tbl2, dist = "beta",
                            hesim_data = hesim_dat),
               row_msg)  
  expect_error(stateval_tbl(tbl2[, !"grp", with = FALSE], dist = "beta",
                            hesim_data = hesim_dat))
  
  # Correct columns for distribution
  ## Beta
  expect_error(stateval_tbl(tbl[,.(strategy_id, grp, state_id, mean)], 
                            dist = "beta", hesim_data = hesim_dat),
               paste0("If a beta distribution is specified, then tbl must ",
               "either contain the columns 'mean' and 'se' or 'shape1' and ", 
               "'shape2'."))
  
  ## Gamma
  expect_error(stateval_tbl(tbl[,.(strategy_id, grp, state_id, mean)], 
                            dist = "gamma", hesim_data = hesim_dat),
               paste0("If a gamma distribution is specified, then tbl must ",
               "either contain the columns 'mean' and 'se', 'shape' and ",
               "'rate', or 'shape' and 'scale'.")) 
  
  ## Lognormal
  expect_error(stateval_tbl(tbl[,.(strategy_id, grp, state_id, mean)], 
                            dist = "lnorm", hesim_data = hesim_dat),
               paste0("If a lognormal distribution is specified, then tbl must ",
               "contain the columns 'meanlog' and 'sdlog'."))   
  
  ## Uniform
  expect_error(stateval_tbl(tbl[,.(strategy_id, grp, state_id, min)], 
                            dist = "unif", hesim_data = hesim_dat),
               paste0("If a uniform distribution is specified, then tbl must ",
                      "contain the columns 'min' and 'max'."))   
})

# StateVals$sim ----------------------------------------------------------------
# Linear model
fit_costs_medical <- stats::lm(costs ~ female + state_name, 
                               data = psm4_exdata$costs$medical)
edat <- expand(hesim_dat, by = c("strategies", "patients", "states"))
costvals_medical <- create_StateVals(fit_costs_medical, input_data = edat, n = N)
costvals_medical$sim(t = c(1, 2, 3), type = "predict")

test_that("StateVals$sim", {
  # Time constant state values 
  expect_equal(c(costvals_medical$input_data$X$mu %*% t(costvals_medical$params$coefs)),
              costvals_medical$sim(t = 2, type = "predict")$value)
  
  ## data must be of class 'input_mats'
  expect_error(StateVals$new(input_mats = 3, params = 2)$sim(t = 2)) 
  
  # Time varying state values
  tbl <- data.table(state_id = rep(c(1, 2, 3), 2),
                    time_start = rep(c(0, 4), each = 3),
                    est = c(1000, 1500, 2000, 2000, 3500, 4000))
  stateval_tbl <- stateval_tbl(tbl, 
                               dist = "fixed", 
                               hesim_data = hesim_dat) 
  mod <- create_StateVals(stateval_tbl, n = 2)
  pred <- mod$sim(t = c(1, 2, 3, 4, 5), type = "predict")
  expect_equal(pred[sample == 1 & strategy_id == 1 & 
                      patient_id == 1 & state_id == 1 & time == 1]$value,
               tbl[state_id == 1 & time_start == 0, est])
  expect_equal(pred[sample == 1 & strategy_id == 1 & 
                      patient_id == 1 & state_id == 3 & time == 5]$value,
               tbl[state_id == 3 & time_start == 4, est])
})

# sim_ev -----------------------------------------------------------------------
n_samples <- 5

# Helper function
R_sim_wlos <- function(prob, stvals, times, dr = .03){
  yvals <- exp(-dr * times) * stvals * prob
  return(pracma::trapz(x = times, y = yvals))
}

wlos_test <- function(econmod, s, k, i, h, dr = .03){
  costs1 <- econmod$costs_[sample == s & strategy_id == k & 
                             patient_id == i & state_id == h]
  stprobs1 <- econmod$stateprobs_[sample == s & strategy_id == k & 
                                  patient_id == i & state_id == h]
  sv_params <- econmod$cost_models[[1]]$params
  sv_obs <- which(sv_params$strategy_id == k & sv_params$patient_id == i &
                    sv_params$state_id == h)
  sv <- sv_params$value[sv_obs, s]
  if (!is.null(sv_params$time_intervals)){
    ti <- findInterval(stprobs1$t, sv_params$time_intervals$time_stop,
                       left.open = TRUE) 
    sv <- rep(sv, table(ti))
  } 
  expect_equal(costs1$costs, 
               R_sim_wlos(stprobs1$prob, sv, stprobs1$t, dr),
               tolerance = .001, scale = 1)
}

# Simulate state probabilities
alpha1 <- rbind(c(200, 400, 600, 800),
                 c(0, 500, 200, 300),
                 c(0, 0, 500, 500),
                 c(0, 0, 0, 1))
alpha2 <- rbind(c(400, 400, 400, 800),
                c(0, 600, 100, 300),
                c(0, 0, 700, 300),
                c(0, 0, 0, 1))
pmat <- data.table(sample = rep(1:n_samples, each = 2),
  strategy_id = rep(1:2, times = 5),
  patient_id = 1,
  rbind(
    rdirichlet_mat(n = 5, alpha = alpha1, output = "data.table"),
    rdirichlet_mat(n = 5, alpha = alpha2, output = "data.table"))
)
transmod <- CohortDtstmTrans$new(params = tparams_transprobs(pmat))
econmod <- CohortDtstm$new(trans_model = transmod)
econmod$sim_stateprobs(n_cycles = 5)

# Simulate costs
hesim_dat <- list(strategies = data.table(strategy_id = 1:2),
                  patients = data.table(patient_id = 1),
                  states = data.table(state_id = 1:3))
stval_tbl <- stateval_tbl(data.table(state_id = 1:3,
                                     mean = c(1701, 1774, 6948),
                                     se = c(1701, 1774, 6948)),
                          dist = "gamma",
                          hesim_data = hesim_dat)
econmod$cost_models <- list(create_StateVals(stval_tbl, n = n_samples))
econmod$sim_costs(dr = .03)

# Tests
wlos_test(econmod, s = 1, k = 1, i = 1, h = 1, dr = .03)
wlos_test(econmod, s = 3, k = 2, i = 1, h = 1, dr = .03)

## Cannot have same discount rate twice
expect_error(econmod$sim_costs(dr = c(.05, .05)))

## Must first simulate stateprobs_
econmod$stateprobs_ <- NULL
expect_error(econmod$sim_costs())

##  Incorrect number of PSA samples
econmod$sim_stateprobs(n_cycles = 5)
econmod$cost_models[[2]] <- create_StateVals(stval_tbl, n = n_samples - 1)
expect_error(econmod$sim_costs())
econmod$cost_models[[2]] <- NULL

## Incorrect number of states
stval_tbl2 <- stateval_tbl(data.table(state_id = 1:4,
                                      est = rep(0, 4)),
                           dist = "fixed",
                           hesim_data = hesim_dat)
econmod$cost_models <- list(create_StateVals(stval_tbl2, n = n_samples))
expect_error(econmod$sim_costs())

## Time varying state values
stval_tbl_tv <- stateval_tbl(data.table(strategy_id = rep(hesim_dat$strategies$strategy_id,
                                                          each = 2),
                                        time_start = c(0, 2, 0, 2),
                                        est = c(1000, 2, 3000, 10)),
                             dist = "fixed",
                             hesim_data = hesim_dat)
econmod$cost_models <- list(create_StateVals(stval_tbl_tv, n = n_samples))
econmod$sim_costs(dr = .03)
wlos_test(econmod, s = 2, k = 2, i = 1, h = 2, dr = .03)

## Using method = "starting" option
stval_tbl_starting <- stateval_tbl(data.table(strategy_id = hesim_dat$strategies$strategy_id,
                                              est = c(1000, 2000)),
                                   dist = "fixed",
                                   hesim_data = hesim_dat)
econmod$cost_models <- list(create_StateVals(stval_tbl_starting, n = n_samples,
                                             method = "starting"))
expect_equal(econmod$cost_models[[1]]$method, "starting")

### With all costs in first health state
econmod$sim_costs()
costs <- dcast(econmod$costs_, sample + strategy_id + patient_id + grp_id +
                 dr + category ~ state_id, value.var = "costs")
expect_true(all(costs[strategy_id == 1][["1"]] == 1000))
expect_true(all(costs[strategy_id == 2][["1"]] == 2000))
expect_true(all(costs[["2"]] == 0))
expect_true(all(costs[["3"]] == 0))

### With costs in 2 health states
econmod$trans_model$start_stateprobs <- c(.5, .5, 0, 0)
econmod$sim_stateprobs(n_cycles = 5)
econmod$sim_costs()
costs <- dcast(econmod$costs_, sample + strategy_id + patient_id + grp_id +
                 dr + category ~ state_id, value.var = "costs")
expect_true(all(costs[strategy_id == 1][["1"]] == 1000 * .5))
expect_true(all(costs[strategy_id == 1][["2"]] == 1000 * .5))
expect_true(all(costs[strategy_id == 2][["1"]] == 2000 * .5))
expect_true(all(costs[strategy_id == 2][["2"]] == 2000 * .5))
expect_true(all(costs[["3"]] == 0))