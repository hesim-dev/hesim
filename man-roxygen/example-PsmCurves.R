library("flexsurv")
N_SAMPLES <- 5 # Number of parameter samples for PSA

# Consider a 3-state model where there is a 
# progression-free survival (PFS) and an
# overall survival (OS) endpoint

# (0) Model setup
hesim_dat <- hesim_data(
  strategies = data.frame(
    strategy_id = c(1, 2),
    strategy_name = c("SOC", "New 1")
  ),
  patients = data.frame(
    patient_id = 1
  )
)

# (1) Parameterize survival models
## (1.1) If patient-level data is available, 
## we can fit survival models

### (1.1.1) Data for estimation (for simplicity, only use 2 strategies)
surv_est_data <- as_pfs_os(
  onc3[strategy_name != "New 2"], 
  patient_vars = c("patient_id", "strategy_name")
)
surv_est_data$strategy_name <- droplevels(surv_est_data$strategy_name)

### (1.1.2) Fit models
fit_pfs <- flexsurvreg(Surv(pfs_time, pfs_status) ~ strategy_name,
                       data = surv_est_data, dist = "exp")
fit_os <- flexsurvreg(Surv(os_time, os_status) ~ strategy_name,
                      data = surv_est_data, dist = "exp")
fits <- flexsurvreg_list(pfs = fit_pfs, os = fit_os)

## (1.2) If patient-level data is NOT available, 
## we can construct the parameter objects "manually"

### (1.2.1) Baseline hazard:
### Assume that we know the (log) rate parameters for both PFS and OS 
### for SOC (i.e., the intercept) and their standard error
logint_pfs_est <- -1.7470900
logint_pfs_se <-  0.03866223
logint_os_est <- -2.7487675
logint_os_se <- 0.04845015

### (1.2.2) Relative treatment effect:
### Assume we know the log hazard ratios (and their standard errors) 
### for comparing the new interventions to the SOC
loghr_pfs_est_new1 <- -0.1772028 
loghr_pfs_se_new1 <- 0.05420119
loghr_os_est_new1 <- -0.1603632
loghr_os_se_new1 <- 0.06948962

### (1.2.3) Create "params_surv_list" object by combining the baseline hazard 
### and relative treatment effects
params <- params_surv_list(
  #### Model for PFS
  pfs = params_surv(
    coefs = list( 
      rate = data.frame( # coefficients predict log rate
        intercept = rnorm(N_SAMPLES, logint_pfs_est, logint_pfs_se),
        new1 = rnorm(N_SAMPLES, loghr_pfs_est_new1, loghr_pfs_se_new1)
      )
    ),
    dist = "exp"
  ),
  
  #### Model for OS
  os = params_surv(
    coefs = list(
      rate = data.frame(
        intercept = rnorm(N_SAMPLES, logint_os_est, logint_os_se),
        new1 = rnorm(N_SAMPLES, loghr_os_est_new1, loghr_os_se_new1)
      )
    ),
    dist = "exp"
  )
)

#### The print (and summary) methods for the "params_surv_list" object will 
#### summarize each of the model terms, which is a good way to check
#### if it's been setup correctly
params 

# (2) Simulation
## (2.1) Construct the model
### (2.1.1) Case where patient-level data was available
### Use create_PsmCurves.params_flexsurvreg_list() method
surv_input_data <- expand(hesim_dat, by = c("strategies", "patients"))
psm_curves1 <- create_PsmCurves(fits, input_data = surv_input_data, 
                                n = N_SAMPLES,
                                uncertainty = "normal",
                                est_data = surv_est_data)

### (2.1.2) Case where patient-level data was NOT available
### Use create_PsmCurves.params_surv_list() method
surv_input_data$intercept <- 1
surv_input_data$new1 <- ifelse(surv_input_data$strategy_name == "New 1", 
                               1, 0)
psm_curves2 <- create_PsmCurves(params, input_data = surv_input_data)

## (2.2) Summarize survival models
## There are minor discrepancies between the case where models were fit
## with flexsurvreg() and the case where the "params_surv_list" object
## was constructed manually due to differences in the random draws
## of the parameter samples. These differences are decreasing in the size 
## of N_SAMPLES
times <- seq(0, 10, 1/12) # Monthly times

### Quantiles
head(psm_curves1$quantile(p = c(.25, .5, .75)))
head(psm_curves2$quantile(p = c(.25, .5, .75)))

### Survival curves
head(psm_curves1$survival(t = times))
head(psm_curves2$survival(t = times))

### Restricted mean survival
head(psm_curves1$rmst(t = c(2, 5)))
head(psm_curves2$rmst(t = c(2, 5)))