library("flexsurv")

# Simulation data
strategies <- data.frame(strategy_id = c(1, 2, 3))
patients <- data.frame(patient_id = seq(1, 3),
                       age = c(45, 50, 60),
                       female = c(0, 0, 1))

# Multi-state model with transition specific models
tmat <- rbind(c(NA, 1, 2),
              c(NA, NA, 3),
              c(NA, NA, NA))
fits <- vector(length = max(tmat, na.rm = TRUE), mode = "list")
for (i in 1:length(fits)){
  fits[[i]] <- flexsurvreg(Surv(years, status) ~ 1,
                           data = bosms3[bosms3$trans == i, ],
                           dist = "exp")
}
fits <- flexsurvreg_list(fits)

# Simulation model
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients)
fits_data <- expand(hesim_dat)
transmod <- create_IndivCtstmTrans(fits, input_data = fits_data,
                                   trans_mat = tmat,
                                   n = 2)
head(transmod$hazard(c(1, 2, 3)))
head(transmod$cumhazard(c(1, 2, 3)))

## Simulate disease progression and state probabilities together
transmod$sim_stateprobs(t = c(0, 5, 10))[t == 5]

## Simulate disease progression and state probabilities separately
disprog <- transmod$sim_disease(max_t = 10)
transmod$sim_stateprobs(t = c(0, 5, 10), disprog = disprog)[t == 5]