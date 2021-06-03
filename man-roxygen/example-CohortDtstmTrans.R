library("msm")
library("data.table")
set.seed(101)

# We consider two examples that have the same treatment strategies and patients.
# One model is parameterized by fitting a multi-state model with the "msm"
# package; in the second model, the parameters are entered "manually" with
# a "params_mlogit_list" object.

# MODEL SETUP
strategies <- data.table(
  strategy_id = c(1, 2, 3),
  strategy_name = c("SOC", "New 1", "New 2")
)
patients <- data.table(patient_id = 1:2)
hesim_dat <- hesim_data(
  strategies = strategies,
  patients = patients
)

# EXAMPLE #1: msm
## Fit multi-state model with panel data via msm
qinit <- rbind(
  c(0, 0.28163, 0.01239),
  c(0, 0, 0.10204),
  c(0, 0, 0)
)
fit <- msm(state_id ~ time, subject = patient_id,
           data = onc3p[patient_id %in% sample(patient_id, 100)],
           covariates = list("1-2" =~ strategy_name),
           qmatrix = qinit)

## Simulation model
transmod_data <- expand(hesim_dat)
transmod <- create_CohortDtstmTrans(fit,
                                    input_data = transmod_data,
                                    cycle_length = 1/2,
                                    fixedpars = 2,
                                    n = 2)
transmod$sim_stateprobs(n_cycles = 2)

# EXAMPLE #2: params_mlogit_list
## Input data
transmod_data[, intercept := 1]
transmod_data[, new1 := ifelse(strategy_name == "New 1", 1, 0)]
transmod_data[, new2 := ifelse(strategy_name == "New 2", 1, 0)]

## Parameters
n <- 10
transmod_params <- params_mlogit_list(
  
  ## Transitions from stable state (stable -> progression, stable -> death)
  stable = params_mlogit(
    coefs = list(
      progression = data.frame(
        intercept = rnorm(n, -0.65, .1),
        new1 = rnorm(n, log(.8), .02),
        new2 = rnorm(n, log(.7, .02))
      ),
      death = data.frame(
        intercept = rnorm(n, -3.75, .1),
        new1 = rep(0, n),
        new2 = rep(0, n)
      )
    )
  ),
  
  ## Transition from progression state (progression -> death)
  progression = params_mlogit(
    coefs = list(
      death = data.frame(
        intercept = rnorm(n, 2.45, .1),
        new1 = rep(0, n),
        new2 = rep(0, n)
      )
    )
  )
)
transmod_params

## Simulation model
tmat <- rbind(c(0, 1, 2),
              c(NA, 0, 1),
              c(NA, NA, NA))
transmod <- create_CohortDtstmTrans(transmod_params, 
                                    input_data = transmod_data,
                                    trans_mat = tmat, cycle_length = 1)
transmod$sim_stateprobs(n_cycles = 2)

\dontshow{
  pb <- expmat(coef(fit)$baseline)[, , 1]
  
  ## From stable
  b1 <- log(pb[1, 2]/(1 - pb[1, 2] - pb[1, 3]))
  b2 <- log(pb[1, 3]/(1 - pb[1, 2] - pb[1, 3]))
  exp(b1)/(1 + exp(b1) + exp(b2))
  exp(b2)/(1 + exp(b1) + exp(b2))
  
  ### From progression
  b <- qlogis(pb[2, 2])
}