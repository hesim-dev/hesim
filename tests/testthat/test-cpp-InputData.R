context("InputData.h unit tests")
library("flexsurv")
library("data.table")

dt.strategies <- data.table(strategy_id = c(1, 2))
dt.patients <- data.table(patient_id = seq(1, 3), 
                          age = c(45, 47, 60),
                          female = c(1, 0, 0),
                          group = factor(c("Good", "Medium", "Poor")))
dt.lines <- lines_dt(list(c(1, 2, 5), c(1, 2)))
dt.states <- data.frame(state_id =  seq(1, 3),
                         state_name = factor(paste0("state", seq(1, 3))))
dt.trans <- data.frame(transition_id = seq(1, 4),
                       from = c(1, 1, 2, 2),
                       to = c(2, 3, 1, 3))
hesim.dat <- hesim_data(strategies = dt.strategies,
                        lines = dt.lines,
                        patients = dt.patients,
                        states = dt.states,
                        transitions = dt.trans)

# InputData --------------------------------------------------------------------
# strategy_id + line + patient_id + state_id
expanded.hesim.dat <- expand_hesim_data(hesim.dat, by = c("strategies", "patients", 
                                                          "lines", "states"))
f <- formula(~ age)
input.dat <- form_input_data(f, expanded.hesim.dat)

test_InputData_row <- function(strategy_id, line, patient_id, state_id){
  obs <- which(input.dat$strategy_id == strategy_id & input.dat$line == line & 
         input.dat$patient_id == patient_id & input.dat$state_id == state_id)
  R.rowvec <- input.dat$X[obs,]
  names(R.rowvec) <- NULL
  C.rowvec <- hesim:::C_test_InputData(input.dat, param_id = 0, strategy_id = strategy_id - 1,
                         line = line - 1, patient_id = patient_id - 1, health_id = state_id - 1)
  expect_equal(R.rowvec, c(C.rowvec))
}

test_InputData_row(strategy_id = 2, line = 1, patient_id = 1, state_id = 3)
test_InputData_row(strategy_id = 2, line = 1, patient_id = 3, state_id = 3)
test_InputData_row(strategy_id = 1, line = 2, patient_id = 2, state_id = 2)
test_InputData_row(strategy_id = 1, line = 3, patient_id = 1, state_id = 1)
test_InputData_row(strategy_id = 1, line = 3, patient_id = 1, state_id = 2)
test_InputData_row(strategy_id = 1, line = 3, patient_id = 2, state_id = 2)

# strategy_id + line + patient_id + transition_id
expanded.hesim.dat <- expand_hesim_data(hesim.dat, by = c("strategies", "patients", 
                                                          "lines", "transitions"))
f <- formula(~ age)
input.dat <- form_input_data(f, expanded.hesim.dat)

test_InputData_row <- function(strategy_id, line, patient_id, transition_id){
  obs <- which(input.dat$strategy_id == strategy_id & input.dat$line == line & 
         input.dat$patient_id == patient_id & input.dat$transition_id == transition_id)
  R.rowvec <- input.dat$X[obs,]
  names(R.rowvec) <- NULL
  C.rowvec <- hesim:::C_test_InputData(input.dat, param_id = 0, strategy_id = strategy_id - 1,
                         line = line - 1, patient_id = patient_id - 1, health_id = transition_id - 1)
  expect_equal(R.rowvec, c(C.rowvec))
}

test_InputData_row(strategy_id = 1, line = 3, patient_id = 1, transition_id = 1)
test_InputData_row(strategy_id = 2, line = 1, patient_id = 1, transition_id = 3)
test_InputData_row(strategy_id = 1, line = 3, patient_id = 2, transition_id = 2)

# strategy_id + line + patient_id
expanded.hesim.dat <- expand_hesim_data(hesim.dat, by = c("strategies", "patients", 
                                                          "lines"))
f <- formula(~ age)
input.dat <- form_input_data(f, expanded.hesim.dat)

test_InputData_row <- function(strategy_id, line, patient_id){
  obs <- which(input.dat$strategy_id == strategy_id & input.dat$line == line & 
         input.dat$patient_id == patient_id)
  R.rowvec <- input.dat$X[obs,]
  names(R.rowvec) <- NULL
  C.rowvec <- hesim:::C_test_InputData(input.dat, param_id = 0, strategy_id = strategy_id - 1,
                         line = line - 1, patient_id = patient_id - 1)
  expect_equal(R.rowvec, c(C.rowvec))
}

test_InputData_row(strategy_id = 1, line = 3, patient_id = 1)
test_InputData_row(strategy_id = 2, line = 1, patient_id = 1)
test_InputData_row(strategy_id = 1, line = 3, patient_id = 2)

# strategy_id + patient_id + state_id
expanded.hesim.dat <- expand_hesim_data(hesim.dat, by = c("strategies", "patients", 
                                                          "states"))
f <- formula(~ age)
input.dat <- form_input_data(f, expanded.hesim.dat)

test_InputData_row <- function(strategy_id, patient_id, state_id){
  obs <- which(input.dat$strategy_id == strategy_id &  input.dat$patient_id == patient_id &
                input.dat$state_id == state_id)
  R.rowvec <- input.dat$X[obs,]
  names(R.rowvec) <- NULL
  C.rowvec <- hesim:::C_test_InputData(input.dat, param_id = 0, strategy_id = strategy_id - 1,
                                        patient_id = patient_id - 1, health_id = state_id - 1)
  expect_equal(R.rowvec, c(C.rowvec))
}

test_InputData_row(strategy_id = 1, patient_id = 1, state_id = 2)
test_InputData_row(strategy_id = 2, patient_id = 1, state_id = 3)
test_InputData_row(strategy_id = 1, patient_id = 2, state_id = 1)

# strategy_id + patient_id 
expanded.hesim.dat <- expand_hesim_data(hesim.dat, by = c("strategies", "patients"))
f <- formula(~ age)
input.dat <- form_input_data(f, expanded.hesim.dat)

test_InputData_row <- function(strategy_id, patient_id){
  obs <- which(input.dat$strategy_id == strategy_id &  input.dat$patient_id == patient_id)
  R.rowvec <- input.dat$X[obs,]
  names(R.rowvec) <- NULL
  C.rowvec <- hesim:::C_test_InputData(input.dat, param_id = 0, strategy_id = strategy_id - 1,
                                        patient_id = patient_id - 1)
  expect_equal(R.rowvec, c(C.rowvec))
}

test_InputData_row(strategy_id = 1, patient_id = 1)
test_InputData_row(strategy_id = 1, patient_id = 2)
test_InputData_row(strategy_id = 1, patient_id = 3)
test_InputData_row(strategy_id = 2, patient_id = 1)
test_InputData_row(strategy_id = 2, patient_id = 2)
test_InputData_row(strategy_id = 2, patient_id = 3)

# TimeFun ----------------------------------------------------------------------
expanded.hesim.dat <- expand_hesim_data(hesim.dat, by = c("strategies", "patients"))
f <- formula(~ age)
input.dat <- form_input_data(f, expanded.hesim.dat)