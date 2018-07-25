context("input_data.h unit tests")
library("flexsurv")
library("data.table")

dt_strategies <- data.table(strategy_id = c(1, 2))
dt_patients <- data.table(patient_id = seq(1, 3), 
                          age = c(45, 47, 60),
                          female = c(1, 0, 0),
                          group = factor(c("Good", "Medium", "Poor")))
dt_lines <- lines_dt(list(c(1, 2, 5), c(1, 2)))
dt_states <- data.frame(state_id =  seq(1, 3),
                         state_name = factor(paste0("state", seq(1, 3))))
dt_trans <- data.frame(transition_id = seq(1, 4),
                       from = c(1, 1, 2, 2),
                       to = c(2, 3, 1, 3))
hesim_dat <- hesim_data(strategies = dt_strategies,
                        lines = dt_lines,
                        patients = dt_patients,
                        states = dt_states,
                        transitions = dt_trans)

# obs_index --------------------------------------------------------------------
## strategy_id + line + patient_id + state_id
expanded_hesim_dat <- expand_hesim_data(hesim_dat, by = c("strategies", "patients", 
                                                          "lines", "states"))
f <- formula_list(mu = ~ age)
input_dat <- form_input_data(f, expanded_hesim_dat)

get_R_rowvec <- function(obs){
  v <- input_dat$X[[1]][obs, ]
  names(v) <- NULL
  return(v)
}

test_obs_index <- function(strategy_id, line, patient_id, state_id){
  R_index <- which(input_dat$strategy_id == strategy_id & input_dat$line == line & 
         input_dat$patient_id == patient_id & input_dat$state_id == state_id)
  C_index <- hesim:::C_test_obs_index(input_dat, strategy_id = strategy_id - 1,
                                     line = line - 1, patient_id = patient_id - 1,
                                     health_id = state_id - 1)
  expect_equal(R_index - 1, C_index)
}

test_obs_index(strategy_id = 2, line = 1, patient_id = 1, state_id = 3)
test_obs_index(strategy_id = 2, line = 1, patient_id = 3, state_id = 3)
test_obs_index(strategy_id = 1, line = 2, patient_id = 2, state_id = 2)
test_obs_index(strategy_id = 1, line = 3, patient_id = 1, state_id = 1)
test_obs_index(strategy_id = 1, line = 3, patient_id = 1, state_id = 2)
test_obs_index(strategy_id = 1, line = 3, patient_id = 2, state_id = 2)

# strategy_id + line + patient_id + transition_id
expanded_hesim_dat <- expand_hesim_data(hesim_dat, by = c("strategies", "patients", 
                                                          "lines", "transitions"))
f <- formula_list(~ age)
input_dat <- form_input_data(f, expanded_hesim_dat)

test_obs_index <- function(strategy_id, line, patient_id, transition_id){
  R_index <- which(input_dat$strategy_id == strategy_id & input_dat$line == line & 
         input_dat$patient_id == patient_id & input_dat$transition_id == transition_id)
  C_index <- hesim:::C_test_obs_index(input_dat, strategy_id = strategy_id - 1,
                         line = line - 1, patient_id = patient_id - 1, health_id = transition_id - 1)
  expect_equal(R_index - 1, C_index)
}

test_obs_index(strategy_id = 1, line = 3, patient_id = 1, transition_id = 1)
test_obs_index(strategy_id = 2, line = 1, patient_id = 1, transition_id = 3)
test_obs_index(strategy_id = 1, line = 3, patient_id = 2, transition_id = 2)

# strategy_id + line + patient_id
expanded_hesim_dat <- expand_hesim_data(hesim_dat, by = c("strategies", "patients", 
                                                          "lines"))
f <- formula_list(~ age)
input_dat <- form_input_data(f, expanded_hesim_dat)

test_obs_index <- function(strategy_id, line, patient_id){
  R_index <- which(input_dat$strategy_id == strategy_id & input_dat$line == line & 
                   input_dat$patient_id == patient_id)
  C_index <- hesim:::C_test_obs_index(input_dat, strategy_id = strategy_id - 1,
                                       line = line - 1, patient_id = patient_id - 1)
  expect_equal(R_index - 1, C_index)
}

test_obs_index(strategy_id = 1, line = 3, patient_id = 1)
test_obs_index(strategy_id = 2, line = 1, patient_id = 1)
test_obs_index(strategy_id = 1, line = 3, patient_id = 2)

# strategy_id + patient_id + state_id
expanded_hesim_dat <- expand_hesim_data(hesim_dat, by = c("strategies", "patients", 
                                                          "states"))
f <- formula_list(~ age)
input_dat <- form_input_data(f, expanded_hesim_dat)

test_obs_index <- function(strategy_id, patient_id, state_id){
  R_index <- which(input_dat$strategy_id == strategy_id &  input_dat$patient_id == patient_id &
                   input_dat$state_id == state_id)
  C_index <- hesim:::C_test_obs_index(input_dat, strategy_id = strategy_id - 1,
                                      patient_id = patient_id - 1, health_id = state_id - 1)
  expect_equal(R_index - 1, C_index)
}

test_obs_index(strategy_id = 1, patient_id = 1, state_id = 2)
test_obs_index(strategy_id = 2, patient_id = 1, state_id = 3)
test_obs_index(strategy_id = 1, patient_id = 2, state_id = 1)

# strategy_id + patient_id 
expanded_hesim_dat <- expand_hesim_data(hesim_dat, by = c("strategies", "patients"))
f <- formula_list(~ age)
input_dat <- form_input_data(f, expanded_hesim_dat)

test_obs_index <- function(strategy_id, patient_id){
  R_index <- which(input_dat$strategy_id == strategy_id &  input_dat$patient_id == patient_id)
  C_index <- hesim:::C_test_obs_index(input_dat, strategy_id = strategy_id - 1,
                                        patient_id = patient_id - 1)
  expect_equal(R_index - 1, C_index)
}

test_obs_index(strategy_id = 1, patient_id = 1)
test_obs_index(strategy_id = 1, patient_id = 2)
test_obs_index(strategy_id = 1, patient_id = 3)
test_obs_index(strategy_id = 2, patient_id = 1)
test_obs_index(strategy_id = 2, patient_id = 2)
test_obs_index(strategy_id = 2, patient_id = 3)

# TimeFun ----------------------------------------------------------------------
# To do