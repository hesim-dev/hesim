library("data.table")

# Input matrices are typically created as part of model objects
# Let's illustrate with a partitioned survival model (PSM)

## Model setup
strategies <- data.frame(strategy_id = c(1, 2),
                         new_strategy = c(0, 1))
patients <- data.frame(patient_id = seq(1, 3),
                       age = c(45, 47, 60),
                       female = c(1, 0, 0),
                       group = factor(c("Good", "Medium", "Poor")))
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients)

## Create survival models for PSM
### Parameters
n <- 2
survmod_params <- params_surv_list(
  # Progression free survival (PFS) 
  pfs = params_surv(
    coefs = list(
      rate = data.frame(intercept = rnorm(n, log(1/5), 1),
                        new_strategy = rnorm(n, log(.8), 1))
      ),
    dist = "exp"
  ),
  
  # Overall survival (OS)
  os = params_surv(
    coefs = list(
      rate = data.frame(intercept = rnorm(n, log(1/10), 1))
    ),
    dist = "exp"
  )
)

### Input data
survmod_input_data <- expand(hesim_dat)[, intercept := 1]

### Model object
survmod <- create_PsmCurves(survmod_params, input_data = survmod_input_data)

## Inspect input data
survmod$input_data # Print "input_mats" object to console
as.data.table(survmod$input_data) # Convert "input_mats" object to data.table