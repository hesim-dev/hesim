context("obs_index.h unit tests")
library("flexsurv")
library("data.table")

strategies <- data.table(strategy_id = c(1, 2))
patients <- data.table(patient_id = seq(1, 3), 
                       age = c(45, 47, 60),
                       female = c(1, 0, 0),
                       group = factor(c("Good", "Medium", "Poor")))
states <- data.frame(state_id =  seq(1, 3),
                     state_name = factor(paste0("state", seq(1, 3))))
trans <- data.frame(transition_id = seq(1, 4),
                    from = c(1, 1, 2, 2),
                    to = c(2, 3, 1, 3))
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients,
                        states = states,
                        transitions = trans)

# obs_index --------------------------------------------------------------------
test_that("obs_index", {

## strategy_id + patient_id + state_id
expanded_hesim_dat <- expand(hesim_dat, by = c("strategies", "patients", 
                                               "states"))
f <- formula_list(mu = ~ age)
X <- create_input_mats(f, expanded_hesim_dat)

get_R_rowvec <- function(obs){
  v <- X$X[[1]][obs, ]
  names(v) <- NULL
  return(v)
}

test_obs_index <- function(strategy_idx, patient_idx, state_idx){
  R_index <- which(X$strategy_id == strategy_idx & X$patient_id == patient_idx & 
                     X$state_id == state_idx)
  C_index <- hesim:::C_test_obs_index(X, strategy_index = strategy_idx - 1,
                                     patient_index = patient_idx - 1,
                                     health_index = state_idx - 1)
  expect_equal(R_index - 1, C_index)
}

test_obs_index(strategy_idx = 2, patient_idx = 1, state_idx = 3)
test_obs_index(strategy_idx = 2, patient_idx = 3, state_idx = 3)
test_obs_index(strategy_idx = 1, patient_idx = 2, state_idx = 2)
test_obs_index(strategy_idx = 1, patient_idx = 1, state_idx = 1)
test_obs_index(strategy_idx = 1, patient_idx = 1, state_idx = 2)
test_obs_index(strategy_idx = 1, patient_idx = 2, state_idx = 2)

# strategy_id + patient_id + transition_id
expanded_hesim_dat <- expand(hesim_dat, by = c("strategies", "patients", 
                                               "transitions"))
f <- formula_list(~ age)
X <- create_input_mats(f, expanded_hesim_dat)

test_obs_index <- function(strategy_idx, patient_idx, transition_idx){
  R_index <- which(X$strategy_id == strategy_idx & X$patient_id == patient_idx & 
                     X$transition_id == transition_idx)
  C_index <- hesim:::C_test_obs_index(X, strategy_index = strategy_idx - 1,
                                      patient_index = patient_idx - 1,
                                      health_index = transition_idx - 1)
  expect_equal(R_index - 1, C_index)
}

test_obs_index(strategy_idx = 1, patient_idx = 1, transition_idx = 1)
test_obs_index(strategy_idx = 2, patient_idx = 1, transition_idx = 3)
test_obs_index(strategy_idx = 1, patient_idx = 2, transition_idx = 2)

# strategy_id + patient_id 
expanded_hesim_dat <- expand(hesim_dat, by = c("strategies", "patients"))
f <- formula_list(~ age)
X <- create_input_mats(f, expanded_hesim_dat)

test_obs_index <- function(strategy_idx, patient_idx){
  R_index <- which(X$strategy_id == strategy_idx &  X$patient_id == patient_idx)
  C_index <- hesim:::C_test_obs_index(X, strategy_index = strategy_idx - 1,
                                        patient_index = patient_idx - 1)
  expect_equal(R_index - 1, C_index)
}

test_obs_index(strategy_idx = 1, patient_idx = 1)
test_obs_index(strategy_idx = 1, patient_idx = 2)
test_obs_index(strategy_idx = 1, patient_idx = 3)
test_obs_index(strategy_idx = 2, patient_idx = 1)
test_obs_index(strategy_idx = 2, patient_idx = 2)
test_obs_index(strategy_idx = 2, patient_idx = 3)

# ID getters  
dat <- expand(hesim_dat, by = c("strategies", "patients", "states"))
X <- input_mats(X = list(mu = model.matrix(~ age, dat)),
                strategy_id = dat$strategy_id,
                n_strategies = length(unique(dat$strategy_id)),
                patient_id = dat$patient_id,
                n_patients = length(unique(dat$patient_id)),
                state_id = dat$state_id,
                n_states = length(unique(dat$state_id)))
strategy_idx <- 1
patient_idx <- 2
health_idx <- 1
obs <- hesim:::C_test_obs_index(X, strategy_index = strategy_idx - 1,
                                patient_index = patient_idx - 1,
                                health_index = health_idx - 1)
test_obs_id <- function(strategy_idx, patient_idx, health_idx, member){
  obs <- hesim:::C_test_obs_index(X, strategy_index = strategy_idx - 1,
                                  patient_index = patient_idx - 1,
                                  health_index = health_idx - 1)
  R_id <- X[[member]][obs + 1]
  C_id <- hesim:::C_test_obs_index_ids(X, strategy_index = strategy_idx - 1,
                                       patient_index = patient_idx - 1,
                                       health_index = health_idx - 1,
                                       member = member)
  expect_equal(R_id, C_id)
}
test_obs_id(strategy_idx = 2, patient_idx = 3, health_idx = 2, "strategy_id")
test_obs_id(strategy_idx = 2, patient_idx = 3, health_idx = 2, "patient_id")
test_obs_id(strategy_idx = 1, patient_idx = 2, health_idx = 2, "patient_id")
test_obs_id(strategy_idx = 1, patient_idx = 2, health_idx = 1, "state_id")
test_obs_id(strategy_idx = 1, patient_idx = 2, health_idx = 2, "state_id")
test_obs_id(strategy_idx = 1, patient_idx = 2, health_idx = 3, "state_id")
}) # end test obs_index