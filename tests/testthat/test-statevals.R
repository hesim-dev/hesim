context("statevals.R unit tests")
library("flexsurv")
library("data.table")
rm(list = ls())

# Simulation
strategies <- data.table(strategy_id = c(1, 2))
patients <- data.table(patient_id = seq(1, 3),
                       grp_id = c(1, 1, 2),
                       age = c(45, 50, 60),
                       female = c(0, 0, 1))
states <- data.frame(state_id =  seq(1, 3),
                     state_name = paste0("state", seq(1, 3)))
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients,
                        states = states)
N <- 5

# State values -----------------------------------------------------------------
# stateval_tbl
test_that("stateval_tbl", {
  n_grps <- 2
  n_strategies <- nrow(strategies)
  n_states <- nrow(states)
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
  expect_equal(nrow(mod$params$mu), nrow(stateval_tbl) * 3 * 2)
  expect_equal(ncol(mod$params$mu), 2)
  expect_true(all(mod$params$mu[c(1, 4, 7, 10), 1] == mod$params$mu[1, 1]))
  
  ## Uniform distribution
  stateval_tbl <- stateval_tbl(tbl[strategy_id == 1 & grp_id == 1,
                                   .(state_id, min, max)], 
                               dist = "unif",  
                               hesim_data = hesim_dat)
  mod <- create_StateVals(stateval_tbl, n = 2)
  expect_true(all(mod$params$mu >= .3 & mod$params$mu <= .5))
  
  ## Normal distribution
  stateval_tbl <- stateval_tbl(tbl[strategy_id == 1 & grp_id == 1,
                                   .(state_id, mean, sd)], 
                               dist = "norm",  
                               hesim_data = hesim_dat)  
  mod <- create_StateVals(stateval_tbl, n = 2) 
  expect_equal(nrow(mod$params$mu), 
               nrow(stateval_tbl) * nrow(patients) * nrow(strategies)) 
  
  # Strategy only
  stateval_tbl <- stateval_tbl(tbl[state_id == 3 & grp_id == 1,
                                   .(strategy_id, mean, se)], 
                               dist = "gamma",
                               hesim_data = hesim_dat)   
  mod <- create_StateVals(stateval_tbl, n = 2)
  expect_equal(nrow(mod$params$mu), nrow(stateval_tbl) * 3 * 3)
  expect_equal(mod$input_mats$state_id, rep(c(1, 2, 3), 3 * 2))
  expect_true(all(mod$params$mu[1:9, 2] == mod$params$mu[1, 2]))
  
  # Strategy, Group, and State
  ## More patients than groups
  stateval_tbl <- stateval_tbl(tbl[, .(strategy_id, grp_id, state_id,
                                       shape1, shape2)], 
                               dist = "beta", 
                               hesim_data = hesim_dat)  
  mod <- create_StateVals(stateval_tbl, n = 1)
  expect_equal(ncol(mod$params$mu), 1)
  expect_equal(mod$input_mats$patient_id,
               rep(rep(patients$patient_id, each = nrow(states)), 
                   nrow(strategies)))
  
  ## Each patient is a unique group
  hesim_dat2 <- hesim_dat
  hesim_dat2$patients <- NULL
  stateval_tbl <- stateval_tbl(tbl, dist = "beta",
                               hesim_data = hesim_dat2)  
  mod <- create_StateVals(stateval_tbl, n = 2)
  expect_equal(nrow(mod$params$mu),
               nrow(stateval_tbl))
  
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
costvals_medical <- create_StateVals(fit_costs_medical, data = edat, n = N)
costvals_medical$sim(t = c(1, 2, 3), type = "predict")

test_that("StateVals$sim", {
  expect_equal(c(costvals_medical$input_mats$X$mu %*% t(costvals_medical$params$coefs)),
              costvals_medical$sim(t = 2, type = "predict")$value)
  
  # data must be of class 'input_mats'
  expect_error(StateVals$new(data = 3, params = 2)$sim(t = 2)) 
})

