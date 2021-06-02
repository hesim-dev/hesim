library("flexsurv")

# Simulation data
strategies <- data.frame(strategy_id = c(1, 2, 3))
patients <- data.frame(patient_id = seq(1, 3),
                       age = c(45, 50, 60),
                       female = c(0, 0, 1))
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients)

# Fit survival models
surv_est_data <- psm4_exdata$survival
fit1 <- flexsurv::flexsurvreg(Surv(endpoint1_time, endpoint1_status) ~ age,
                              data = surv_est_data, dist = "exp")
fit2 <- flexsurv::flexsurvreg(Surv(endpoint2_time, endpoint2_status) ~ age,
                              data = surv_est_data, dist = "exp")
fit3 <- flexsurv::flexsurvreg(Surv(endpoint3_time, endpoint3_status) ~ age,
                              data = surv_est_data, dist = "exp")
fits <- flexsurvreg_list(fit1, fit2, fit3)

# Form PsmCurves
surv_input_data <- expand(hesim_dat, by = c("strategies", "patients"))
psm_curves <- create_PsmCurves(fits, input_data = surv_input_data, n = 3,
                               uncertainty = "bootstrap", est_data = surv_est_data)

# Summarize survival curves
head(psm_curves$quantile(p = c(.25, .5, .75)))
head(psm_curves$survival(t = seq(0, 3, by = .1)))
head(psm_curves$rmst(t = c(2, 5)))
