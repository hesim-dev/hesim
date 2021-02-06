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

# stateval_tbl and create_StateVals.stateval_tbl -------------------------------
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

test_that(paste0("create_StateVals.stateval_tbl() returns correct output with state_id ",
                 "column and beta distribution"), {
  sv <- stateval_tbl(
    tbl[strategy_id == 1 & grp == 1, .(state_id, mean, se)], 
    dist = "beta"
  )
  mod <- create_StateVals(sv, n = 2, hesim_data = hesim_dat)
  expect_equal(nrow(mod$params$value), nrow(sv) * 3 * 2)
  expect_equal(ncol(mod$params$value), 2)
  expect_true(all(mod$params$value[c(1, 4, 7, 10), 1] == mod$params$value[1, 1]))
})

test_that("'time_reset' argument is passed correctly with create_StateVals.stateval_tbl()", {
  sv <- stateval_tbl(
    tbl[strategy_id == 1 & grp == 1, .(state_id, mean, se)], 
    dist = "beta"
  )
  mod <- create_StateVals(sv, n = 2, hesim_data = hesim_dat, time_reset = TRUE)
  expect_true(mod$time_reset)
}) 

test_that(paste0("create_StateVals.stateval_tbl() returns correct output with state_id ",
                 "column and uniform distribution"), {
  sv <- stateval_tbl(
    tbl[strategy_id == 1 & grp == 1, .(state_id, min, max)], 
    dist = "unif"
  )
  mod <- create_StateVals(sv, n = 2, hesim_data = hesim_dat)
  expect_true(all(mod$params$value >= .3 & mod$params$value <= .5))  
})

test_that(paste0("create_StateVals.stateval_tbl() returns correct output with state_id ",
                 "column and normal distribution"), {
  sv <- stateval_tbl(
    tbl[strategy_id == 1 & grp == 1, .(state_id, mean, sd)], 
    dist = "norm"
  )  
  mod <- create_StateVals(sv, n = 2, hesim_data = hesim_dat) 
  expect_equal(nrow(mod$params$value), 
               nrow(sv) * n_patients * n_strategies)  
})

test_that(paste0("create_StateVals.stateval_tbl() works as expected with state_id", 
                 " column and custom distribution"), {
  tbl2 <- data.table(state_id = rep(states$state_id, 2),
                    sample = rep(c(1, 2), each = n_states),
                    value = rnorm(6, 4))
  sv <- stateval_tbl(tbl2, dist = "custom")
  
  # Should be error if the distribution is not custom
  expect_error(
    stateval_tbl(tbl2),
    "If 'sample' is in 'tbl', then 'dist' must equal 'custom'"
  )
  
  # n is less than the number of samples in tbl
  expect_warning(mod <- create_StateVals(sv, n = 3, hesim_data = hesim_dat))
  expect_equal(ncol(mod$params$value), 3)
  expect_true(mod$params$value[1, 1] == tbl2[state_id == 1 & sample == 1, value] |
                mod$params$value[1, 1] == tbl2[state_id == 1 & sample == 2, value]) 
  
  # n equals number of samples in tbl
  mod <- create_StateVals(sv, n = 2, hesim_data = hesim_dat) 
  expect_equal(ncol(mod$params$value), 2)
  expect_equal(nrow(mod$params$value), 
               n_states * n_patients * n_strategies) 
  expect_equal(mod$params$value[5, 1],
               tbl2[state_id == 2 & sample == 1, value])
  
  # n is less than the number of samples in tbl
  mod <- create_StateVals(sv, n = 1, hesim_data = hesim_dat) 
  expect_equal(ncol(mod$params$value), 1)
  expect_true(mod$params$value[3, 1] == tbl2[state_id == 3 & sample == 1, value] |
                mod$params$value[3, 1] == tbl2[state_id == 3 & sample == 2, value]) 
})

test_that("create_StateVals.stateval_tbl() returns correct output state_id column and gamma distribution", {
  sv <- stateval_tbl(
    tbl[state_id == 3 & grp == 1, .(strategy_id, mean, se)], 
    dist = "gamma"
  )
  mod <- create_StateVals(sv, n = 2, hesim_data = hesim_dat)
  expect_equal(nrow(mod$params$value), nrow(sv) * 3 * 3)
  expect_equal(mod$params$state_id, rep(c(1, 2, 3), 3 * 2))
  expect_true(all(mod$params$value[1:9, 2] == mod$params$value[1, 2]))
})

test_that(paste0("create_StateVals.stateval_tbl() returns correct output with" ,
                 "patient_id column and gamma distribution"), {
  tbl[, patient_id := grp]
  sv <- stateval_tbl(
    tbl[state_id == 3 & grp == 1, .(strategy_id, mean, se)], 
    dist = "gamma"
  )
  mod <- create_StateVals(sv, n = 2, hesim_data = hesim_dat)
  expect_equal(nrow(mod$params$value), nrow(sv) * 3 * 3)
  expect_equal(mod$params$state_id, rep(c(1, 2, 3), 3 * 2))
  expect_true(all(mod$params$value[1:9, 2] == mod$params$value[1, 2]))
})

test_that("create_StateVals.stateval_tbl() returns correct output when using grp_var", {
  sv <- stateval_tbl(
    tbl[, .(strategy_id, grp, state_id, shape1, shape2)], 
    dist = "beta", 
    grp_var = "grp"
  )  
  mod <- create_StateVals(sv, n = 1, hesim_data = hesim_dat)
  expect_equal(ncol(mod$params$value), 1)
  expect_equal(mod$params$patient_id,
               rep(rep(patients$patient_id, each = nrow(states)), 
                   nrow(strategies))) 
})
  
test_that("create_StateVals.stateval_tbl() returns correct output with state_id and time columns", {
  tbl2 <- tbl[strategy_id == 1]
  tbl2[, time_start := ifelse(grp == 1, 0, 4)]
  tbl2[, time_stop := ifelse(grp == 1, 4, 10)]
  sv <- stateval_tbl(
    tbl2[,.(state_id, time_start, time_stop, mean, se)], 
    dist = "beta"
  )
  mod <- create_StateVals(sv, n = 2, hesim_data = hesim_dat)
  expect_equal(nrow(mod$params$value),
               n_states * n_strategies * n_patients * 2)
})

test_that(paste0("create_StateVals.stateval_tbl() returns an error if the 'states'",
                  "table is not included in hesim_data when it should be"), {
  sv <- stateval_tbl(
    tbl[state_id == 1 & grp == 1, .(strategy_id, mean, se)], 
    dist = "beta"
  )
  expect_error(
    create_StateVals(sv, n = 2),
    paste0("If 'state_id' is not a column in 'object', then ",
           "'hesim_data' must be included as an argument and ", 
           "'states' must be an element of 'hesim_data'.")
    ) 
})

test_that(paste0("create_StateVals.stateval_tbl() returns an error if the 'patients'",
                 "table is not included in hesim_data when it should be"), {
  sv <- stateval_tbl(
   tbl[patient_id == 1 & grp == 1, .(strategy_id, state_id, mean, se)], 
   dist = "beta"
  )
  expect_error(
   create_StateVals(sv, n = 2),
   paste0("If 'patient_id' is not a column in 'object', then ",
          "'hesim_data' must be included as an argument and ", 
          "'patients' must be an element of 'hesim_data'.")
  ) 
})

test_that("stateval_tbl() returns error if the number of rows is wreong", {
  tbl2 <- tbl[,.(strategy_id, state_id, grp, mean, se)]
  tbl2 <- rbind(tbl2, data.table(strategy_id = 2, state_id = 1, grp = 1, 
                                 mean = .5, se = .1))
  expect_error(
    stateval_tbl(tbl2, dist = "beta"),
    paste0("There must only be one row for each combination ",
           "of strategy_id and state_id in 'tbl'.")
  )
})


test_that("stateval_tbl() returns error if the columns aren't correct for the distribution", {
  
  # Beta
  expect_error(stateval_tbl(tbl[,.(strategy_id, grp, state_id, mean)], 
                            dist = "beta"),
               paste0("If a beta distribution is specified, then tbl must ",
               "either contain the columns 'mean' and 'se' or 'shape1' and ", 
               "'shape2'."))
  
  # Gamma
  expect_error(stateval_tbl(tbl[,.(strategy_id, grp, state_id, mean)], 
                            dist = "gamma"),
               paste0("If a gamma distribution is specified, then tbl must ",
               "either contain the columns 'mean' and 'se', 'shape' and ",
               "'rate', or 'shape' and 'scale'.")) 
  
  # Lognormal
  expect_error(stateval_tbl(tbl[,.(strategy_id, grp, state_id, mean)], 
                            dist = "lnorm"),
               paste0("If a lognormal distribution is specified, then tbl must ",
               "contain the columns 'meanlog' and 'sdlog'."))   
  
  # Uniform
  expect_error(stateval_tbl(tbl[,.(strategy_id, grp, state_id, min)], 
                            dist = "unif"),
               paste0("If a uniform distribution is specified, then tbl must ",
                      "contain the columns 'min' and 'max'."))   
})

test_that("stateval_tbl() returns error if the number of rows is incorrect", {
  expect_error(
    stateval_tbl(
      tbl[strategy_id == 1][!(grp ==2 & state_id == 3)],
      dist = "beta"
    ),
    paste0("The number of rows in 'object' should equal 6 which is the number of",
            " unique values of strategy_id, state_id, and patient_id in 'object'.")
  )
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
  sv <- stateval_tbl(tbl, dist = "fixed") 
  mod <- create_StateVals(sv, n = 2, hesim_data = hesim_dat)
  pred <- mod$sim(t = c(1, 2, 3, 4, 5), type = "predict")
  expect_equal(pred[sample == 1 & strategy_id == 1 & 
                      patient_id == 1 & state_id == 1 & time == 1]$value,
               tbl[state_id == 1 & time_start == 0, est])
  expect_equal(pred[sample == 1 & strategy_id == 1 & 
                      patient_id == 1 & state_id == 3 & time == 5]$value,
               tbl[state_id == 3 & time_start == 4, est])
})

# sim_costs --------------------------------------------------------------------
# Helper functions
R_sim_wlos <- function(prob, stvals, times, dr = .03){
  yvals <- exp(-dr * times) * stvals * prob
  return(pracma::trapz(x = times, y = yvals))
}

test_wlos <- function(econmod, s, k, i, h, dr = .03){
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

# Construct model + simulate outcomes
n_samples <- 5
hesim_dat <- list(
  strategies = data.table(strategy_id = 1:2),
  patients = data.table(patient_id = 1),
  states = data.table(state_id = 1:3)
)

## Transition model
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

## Cost model
cost_tbl <- stateval_tbl(data.table(state_id = 1:3,
                                     mean = c(1701, 1774, 6948),
                                     se = c(1701, 1774, 6948)),
                          dist = "gamma")
costmods <- list(create_StateVals(cost_tbl, n = n_samples,
                                  hesim_data = hesim_dat))

## Full economic model + simulation of state probabilities
econmod <- CohortDtstm$new(trans_model = transmod,
                           cost_models = costmods)
econmod$sim_stateprobs(n_cycles = 5)

# Run tests
test_that("sim_costs() produces correct result", {
  econmod$sim_costs(dr = .03)
  test_wlos(econmod, s = 1, k = 1, i = 1, h = 1, dr = .03)
  test_wlos(econmod, s = 3, k = 2, i = 1, h = 1, dr = .03)
})

test_that("sim_costs() cannot have same discount rate twice", {
  expect_error(econmod$sim_costs(dr = c(.05, .05)))
})

test_that("sim_costs() throws error if stateprobs_ have not been simulated", {
  econmod$stateprobs_ <- NULL
  expect_error(econmod$sim_costs())
})

test_that("sim_costs() throws error with incorrect number of PSA samples", {
  econmod$sim_stateprobs(n_cycles = 5)
  econmod$cost_models[[2]] <- create_StateVals(cost_tbl, n = n_samples - 1,
                                               hesim_data = hesim_dat)
  expect_error(
    econmod$sim_costs(),
    paste0("The number of samples in each 'StateVals' model must equal the ",
           "number of samples in the 'stateprobs' object, which is 5.")
  )
})

test_that("sim_costs() throws error with incorrect number of strategies", {
  econmod$cost_models[[2]] <- costmods[[1]]$clone()
  econmod$cost_models[[2]]$params$n_strategies <- 1
  expect_error(
    econmod$sim_costs(),
    paste0("The number of strategies in each 'StateVals' model must equal the ",
           "number of strategies in the 'stateprobs' object, which is 2.")
  )
})

test_that("sim_costs() throws error with incorrect number of patients", {
  econmod$cost_models[[2]] <- costmods[[1]]$clone()
  econmod$cost_models[[2]]$params$n_patients <- 2
  expect_error(
    econmod$sim_costs(),
    paste0("The number of patients in each 'StateVals' model must equal the ",
           "number of patients in the 'stateprobs' object, which is 1.")
  )
})

test_that("sim_costs() throws error with incorrect number of health states", {
  cost_tbl2 <- stateval_tbl(
    data.table(state_id = 1:4,
                est = rep(0, 4)),
                dist = "fixed")
  econmod$cost_models <- list(create_StateVals(cost_tbl2, n = n_samples, 
                                               hesim_data = hesim_dat))
  expect_error(
    econmod$sim_costs(),
    paste0("The number of states in each 'StateVals' model must be one less ",
           "(since state values are not applied to the death state) than ",
           "the number of states in the 'stateprobs' object, which is 4."),
    fixed = TRUE
    )
})

test_that("sim_costs produces correct result with time-varying state values", {
  cost_tbl_tv <- stateval_tbl(
    data.table(strategy_id = rep(hesim_dat$strategies$strategy_id,each = 2),
              time_start = c(0, 2, 0, 2),
              est = c(1000, 2, 3000, 10)),
              dist = "fixed"
  )
  econmod$cost_models <- list(create_StateVals(cost_tbl_tv, n = n_samples,
                                               hesim_data = hesim_dat))
  econmod$sim_costs(dr = .03)
  test_wlos(econmod, s = 2, k = 2, i = 1, h = 2, dr = .03)
})

cost_tbl_starting <- stateval_tbl(
  data.table(strategy_id = hesim_dat$strategies$strategy_id,
  est = c(1000, 2000)),
  dist = "fixed")
econmod$cost_models <- list(
  create_StateVals(cost_tbl_starting, n = n_samples, method = "starting",
                   hesim_data = hesim_dat)
)

test_that("create_StateVals() correctly passes method = 'starting'", {
  expect_equal(econmod$cost_models[[1]]$method, "starting")
})

test_that("sim_costs produces correct result with method = 'starting' and 
          all costs in the first health state", {
  econmod$sim_costs()
  costs <- dcast(econmod$costs_, sample + strategy_id + patient_id + grp_id +
                   dr + category ~ state_id, value.var = "costs")
  expect_true(all(costs[strategy_id == 1][["1"]] == 1000))
  expect_true(all(costs[strategy_id == 2][["1"]] == 2000))
  expect_true(all(costs[["2"]] == 0))
  expect_true(all(costs[["3"]] == 0))
})

test_that("sim_costs produces correct result with method = 'starting' and 
          costs in 2 health states", {
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
})

# sim_qalys --------------------------------------------------------------------
# Add utility model to economic model above
utility_tbl <- stateval_tbl(
  data.table(state_id = 1:3,
              mean = c(1, .8, .7),
              se = c(0, .02, .03)),
              dist = "beta")
utilitymod <- create_StateVals(utility_tbl, n = n_samples, hesim_data = hesim_dat)
econmod <- CohortDtstm$new(trans_model = transmod,
                           utility_model = utilitymod)

test_that("sim_qalys() returns an error if $_stateprobs is NULL", {
  expect_error(
    econmod$sim_qalys(),
    "You must first simulate state probabilities using '$sim_stateprobs'.",
    fixed = TRUE
  )
})

econmod$sim_stateprobs(n_cycles = 5)

test_that("sim_qalys does not return lys column if lys = FALSE", {
  econmod$sim_qalys(lys = FALSE)
  expect_true(!"lys" %in% colnames(econmod$qalys_))
})

test_that("sim_qalys produces correct value of lys", {
  dr <- .03
  econmod$sim_qalys(lys = TRUE, dr = dr, integrate_method = "riemann_left")
  max_t <- max(econmod$stateprobs_$t)
  max_state_id <- max(econmod$stateprobs_$state_id)
  lys <- econmod$stateprobs_[t != max_t & state_id != max_state_id , 
                      .(lys = sum(exp(t * -dr) * prob)), 
                      by = c("sample", "strategy_id", "patient_id", 
                             "grp_id", "state_id")]
  expect_equal(econmod$qalys_$lys, lys$lys)
})