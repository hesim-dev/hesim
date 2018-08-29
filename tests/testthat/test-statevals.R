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
# stateval_means object
test_that("stateval_means", {
  gamma_params <- mom_gamma(c(5000, 7000), c(1000, 1200))
  n <- 6
  vals <- matrix(rgamma(2 * n, 
                        shape = gamma_params$shape, 
                        scale = gamma_params$scale),
                 nrow = n, ncol = 2, byrow = TRUE)
  stval_ests <- stateval_means(values = vals,
                             strategy_id = c(1, 2),
                             patient_id = c(1, 2, 3))
  statevals <- create_StateVals(stval_ests)
  statevals$check()
  predict <- statevals$sim(t = c(1, 2), type = "predict")
  random <- statevals$sim(t = c(1, 2), type = "random")
  expect_equal(predict, random)
  expect_equal(predict[sample == 1 & time == 1]$value,
               statevals$params$mu[, 1])
})

# Linear model
fit_costs_medical <- stats::lm(costs ~ female + state_name, 
                               data = psm4_exdata$costs$medical)
edat <- expand_hesim_data(hesim_dat, by = c("strategies", "patients", "states"))
costvals_medical <- create_StateVals(fit_costs_medical, data = edat, n = N)
costvals_medical$sim(t = c(1, 2, 3), type = "predict")

test_that("StateVals$sim", {
  expect_equal(c(costvals_medical$data$X$mu %*% t(costvals_medical$params$coefs)),
              costvals_medical$sim(t = 2, type = "predict")$value)
  
  # data must be of class 'input_data'
  expect_error(StateVals$new(data = 3, params = 2)$sim(t = 2)) 
  
  # class of 'params' is not supported
  input_dat <- create_input_data(formula_list(~1), edat)
  expect_error(StateVals$new(data = input_dat, params = 2)$sim(3))
})

