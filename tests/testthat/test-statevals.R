context("statevals.R unit tests")
library("flexsurv")
library("data.table")

# Simulation
dt_strategies <- data.table(strategy_id = c(1, 2, 3))
dt_patients <- data.table(patient_id = seq(1, 3),
                          age = c(45, 50, 60),
                          female = c(0, 0, 1))
dt_states <- data.frame(state_id =  seq(1, 3),
                        state_name = paste0("state", seq(1, 3)))
hesim_dat <- hesim_data(strategies = dt_strategies,
                        patients = dt_patients,
                        states = dt_states)
N <- 5

# State values -----------------------------------------------------------------
fit_costs_medical <- stats::lm(costs ~ female + state_name, 
                               data = psm4_exdata$costs$medical)
edat <- expand_hesim_data(hesim_dat, by = c("strategies", "patients", "states"))
costvals_medical <- form_StateVals(fit_costs_medical, data = edat, n = N)
costvals_medical$sim(t = c(1, 2, 3), type = "predict")

test_that("StateVals$sim", {
  expect_equal(c(costvals_medical$data$X$mu %*% t(costvals_medical$params$coefs)),
              costvals_medical$sim(t = 2, type = "predict")$value)
  
  # data must be of class 'input_data'
  expect_error(StateVals$new(data = 3, params = 2)$sim(t = 2)) 
  
  # class of 'params' is not supported
  input_dat <- form_input_data(formula_list(~1), edat)
  expect_error(StateVals$new(data = input_dat, params = 2)$sim(3))
})

