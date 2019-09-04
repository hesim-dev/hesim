context("statevals.R unit tests")
library("data.table")
rm(list = ls())

# Simulation
strategies <- data.table(strategy_id = c(1, 2))
n_strategies <- nrow(strategies)
patients <- data.table(patient_id = seq(1, 3),
                       grp_id = c(1, 1, 2),
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

# State values -----------------------------------------------------------------
# stateval_tbl
test_that("stateval_tbl", {
  n_grps <- 2
  strategy_id <- rep(strategies$strategy_id, each = n_grps * n_states)
  state_id <- rep(states$state_id, times = n_grps * n_strategies)
  grp_id <- rep(rep(c(1, 2), each = n_states), times = n_strategies)
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
  tbl <- data.table(strategy_id = strategy_id, grp_id = grp_id,
                    state_id = state_id, mean = mean, se = se,
                    sd = se,
                    shape1 = shape1,
                    shape2 = shape2, min = .3, max = .5)
  
  # State only
  ## Beta distribution
  stateval_tbl <- stateval_tbl(tbl[strategy_id == 1 & grp_id == 1,
                                   .(state_id, mean, se)], 
                               dist = "beta", 
                               hesim_data = hesim_dat) 
  mod <- create_StateVals(stateval_tbl, n = 2)
  expect_equal(nrow(mod$params$value), nrow(stateval_tbl) * 3 * 2)
  expect_equal(ncol(mod$params$value), 2)
  expect_true(all(mod$params$value[c(1, 4, 7, 10), 1] == mod$params$value[1, 1]))
  
  ### time_reset = TRUE
  mod <- create_StateVals(stateval_tbl, n = 2, time_reset = TRUE)
  expect_true(mod$time_reset)
  
  ## Uniform distribution
  stateval_tbl <- stateval_tbl(tbl[strategy_id == 1 & grp_id == 1,
                                   .(state_id, min, max)], 
                               dist = "unif",  
                               hesim_data = hesim_dat)
  mod <- create_StateVals(stateval_tbl, n = 2)
  expect_true(all(mod$params$value >= .3 & mod$params$value <= .5))
  
  ## Normal distribution
  stateval_tbl <- stateval_tbl(tbl[strategy_id == 1 & grp_id == 1,
                                   .(state_id, mean, sd)], 
                               dist = "norm",  
                               hesim_data = hesim_dat)  
  mod <- create_StateVals(stateval_tbl, n = 2) 
  expect_equal(nrow(mod$params$value), 
               nrow(stateval_tbl) * n_patients * n_strategies) 
  
  ## Custom distribution
  tbl2 <- data.table(state_id = rep(states$state_id, 2),
                    sample = rep(c(1, 2), each = n_states),
                    value = rnorm(6, 4))
  expect_error(stateval_tbl(tbl2, 
                            hesim_data = hesim_dat))
  stateval_tbl <- stateval_tbl(tbl2, 
                               dist = "custom",
                               hesim_data = hesim_dat)
  
  ### n is greater than number of samples in tbl
  expect_warning(mod <- create_StateVals(stateval_tbl, n = 3))
  expect_equal(ncol(mod$params$value), 3)
  expect_true(mod$params$value[1, 1] == tbl2[state_id == 1 & sample == 1, value] |
              mod$params$value[1, 1] == tbl2[state_id == 1 & sample == 2, value]) 
  
  ### n equals number of samples in tbl
  mod <- create_StateVals(stateval_tbl, n = 2) 
  expect_equal(ncol(mod$params$value), 2)
  expect_equal(nrow(mod$params$value), 
               n_states * n_patients * n_strategies) 
  expect_equal(mod$params$value[5, 1],
               tbl2[state_id == 2 & sample == 1, value])
  
  ### n is less than the number of samples in tbl
  mod <- create_StateVals(stateval_tbl, n = 1) 
  expect_equal(ncol(mod$params$value), 1)
  expect_true(mod$params$value[3, 1] == tbl2[state_id == 3 & sample == 1, value] |
              mod$params$value[3, 1] == tbl2[state_id == 3 & sample == 2, value]) 
  
  # Strategy only
  stateval_tbl <- stateval_tbl(tbl[state_id == 3 & grp_id == 1,
                                   .(strategy_id, mean, se)], 
                               dist = "gamma",
                               hesim_data = hesim_dat)   
  mod <- create_StateVals(stateval_tbl, n = 2)
  expect_equal(nrow(mod$params$value), nrow(stateval_tbl) * 3 * 3)
  expect_equal(mod$params$state_id, rep(c(1, 2, 3), 3 * 2))
  expect_true(all(mod$params$value[1:9, 2] == mod$params$value[1, 2]))
  
  # Strategy, Group, and State
  ## More patients than groups
  stateval_tbl <- stateval_tbl(tbl[, .(strategy_id, grp_id, state_id,
                                       shape1, shape2)], 
                               dist = "beta", 
                               hesim_data = hesim_dat)  
  mod <- create_StateVals(stateval_tbl, n = 1)
  expect_equal(ncol(mod$params$value), 1)
  expect_equal(mod$params$patient_id,
               rep(rep(patients$patient_id, each = nrow(states)), 
                   nrow(strategies)))
  
  ## Each patient is a unique group
  hesim_dat2 <- hesim_dat
  hesim_dat2$patients <- NULL
  stateval_tbl <- stateval_tbl(tbl, dist = "beta",
                               hesim_data = hesim_dat2)  
  mod <- create_StateVals(stateval_tbl, n = 2)
  expect_equal(nrow(mod$params$value),
               nrow(stateval_tbl))
  
  # State and time period
  tbl2 <- tbl[strategy_id == 1]
  tbl2[, time_start := ifelse(grp_id == 1, 0, 4)]
  tbl2[, time_stop := ifelse(grp_id == 1, 4, 10)]
  stateval_tbl <- stateval_tbl(tbl2[,.(state_id, time_start, time_stop,
                                      mean, se)], 
                                dist = "beta",
                               hesim_data = hesim_dat)
  mod <- create_StateVals(stateval_tbl, n = 2)
  expect_equal(nrow(mod$params$value),
               n_states * n_strategies * n_patients * 2)
  
  # Errors
  ## state_id and strategy_id are correctly specified
  expect_error(stateval_tbl(tbl[,.(strategy_id, grp_id, mean)], 
                      dist = "beta"))
  expect_error(stateval_tbl(tbl[,.(state_id, grp_id, mean)], 
                      dist = "beta"))  
  expect_error(stateval_tbl(tbl[,.(state_id, strategy_id, grp_id, mean)], 
                      dist = "beta"))  
  
  ## Correct number of rows
  tbl2 <- tbl[,.(strategy_id, state_id, grp_id, mean, se)]
  tbl2 <- rbind(tbl2, data.table(strategy_id = 2, state_id = 1, grp_id = 1, 
                                 mean = .5, se = .1))
  expect_error(stateval_tbl(tbl2, dist = "beta",
                            hesim_data = hesim_dat))
  
  tbl2 <- tbl2[(strategy_id == 1 & state_id == 1) |
               (strategy_id == 2 & state_id == 2)]
  expect_error(stateval_tbl(tbl2, dist = "beta",
                            hesim_data = hesim_dat))  
  expect_error(stateval_tbl(tbl2[, !"grp_id", with = FALSE], dist = "beta",
                            hesim_data = hesim_dat))
  
  ## Correct columns for distribution
  ### Beta
  expect_error(stateval_tbl(tbl[,.(strategy_id, grp_id, state_id, mean)], 
                                dist = "beta", hesim_data = hesim_dat))
  
  ### Gamma
  expect_error(stateval_tbl(tbl[,.(strategy_id, grp_id, state_id, mean)], 
                                dist = "gamma", hesim_data = hesim_dat)) 
  
  ### Lognormal
  expect_error(stateval_tbl(tbl[,.(strategy_id, grp_id, state_id, mean)], 
                                dist = "lnorm", hesim_data = hesim_dat))   
  
  ### Uniform
  expect_error(stateval_tbl(tbl[,.(strategy_id, grp_id, state_id, min)], 
                                dist = "unif", hesim_data = hesim_dat))   
})

# Linear model
fit_costs_medical <- stats::lm(costs ~ female + state_name, 
                               data = psm4_exdata$costs$medical)
edat <- expand(hesim_dat, by = c("strategies", "patients", "states"))
costvals_medical <- create_StateVals(fit_costs_medical, input_data = edat, n = N)
costvals_medical$sim(t = c(1, 2, 3), type = "predict")

test_that("StateVals$sim", {
  # Time constant state values 
  expect_equal(c(costvals_medical$input_mats$X$mu %*% t(costvals_medical$params$coefs)),
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

