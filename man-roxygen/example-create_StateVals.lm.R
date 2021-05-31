library("data.table")
set.seed(10)

# EXAMPLE FOR `create_statevals.lm()`
## Simple example comparing two treatment strategies where
## medical costs vary by sex and health state

## Setup model
hesim_dat <- hesim_data(
  strategies = data.frame(strategy_id = c(1, 2)),
  patients = data.frame(
    patient_id = c(1, 2),
    female = c(1, 0)
  ),
  states = data.frame(
    state_id = c(1, 2, 3),
    state_name = c("state1", "state2", "state3")
  )
)

## Fit model
medcost_estimation_data <- data.table(psm4_exdata$costs$medical)
medcost_estimation_data[, time5 := rbinom(.N, 1, .5)] # Illustrative time dummy
medcost_fit <- lm(costs ~ female + state_name + time5, 
                  data = medcost_estimation_data)

## Create medical cost model
### Allow medical costs to vary across time in addition to by patient and 
### health state
medcost_times <- time_intervals(
  data.frame(time_start = c(0, 3, 5),
            time5 = c(0, 0, 1)) # Time dummy corresponds to time > 5
)
medcost_input_data <- expand(hesim_dat, 
                             by = c("strategies", "patients", "states"),
                             times = medcost_times)
medcost_model <- create_StateVals(medcost_fit, medcost_input_data,
                                  n = 1)

## Explore predictions from medical cost model
### We can assess predictions at multiple time points
medcost_model$sim(t = c(1, 6), type = "predict")