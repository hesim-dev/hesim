library("flexsurv")
library("ggplot2")
theme_set(theme_bw())

# Model setup
strategies <- data.frame(strategy_id = c(1, 2, 3),
                         strategy_name = paste0("Strategy ", 1:3))
patients <- data.frame(patient_id = seq(1, 3),
                       age = c(45, 50, 60),
                       female = c(0, 0, 1))
states <- data.frame(state_id =  seq(1, 3),
                     state_name = paste0("State ", seq(1, 3)))
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients,
                        states = states)
labs <- c(
  get_labels(hesim_dat),
  list(curve = c("Endpoint 1" = 1,
                 "Endpoint 2" = 2,
                 "Endpoint 3" = 3))
)
n_samples <- 2

# Survival models
surv_est_data <- psm4_exdata$survival
fit1 <- flexsurvreg(Surv(endpoint1_time, endpoint1_status) ~ factor(strategy_id),
                    data = surv_est_data, dist = "exp")
fit2 <- flexsurvreg(Surv(endpoint2_time, endpoint2_status) ~ factor(strategy_id),
                    data = surv_est_data, dist = "exp")
fit3 <- flexsurvreg(Surv(endpoint3_time, endpoint3_status) ~ factor(strategy_id),
                    data = surv_est_data, dist = "exp")
fits <- flexsurvreg_list(fit1, fit2, fit3)

surv_input_data <- expand(hesim_dat, by = c("strategies", "patients"))
psm_curves <- create_PsmCurves(fits, input_data = surv_input_data,
                               uncertainty = "bootstrap", est_data = surv_est_data,
                               n = n_samples)

# Cost model(s)
cost_input_data <- expand(hesim_dat, by = c("strategies", "patients", "states"))
fit_costs_medical <- lm(costs ~ female + state_name,
                        data = psm4_exdata$costs$medical)
psm_costs_medical <- create_StateVals(fit_costs_medical,
                                      input_data = cost_input_data,
                                      n = n_samples)

# Utility model
utility_tbl <- stateval_tbl(tbl = data.frame(state_id = states$state_id,
                                             min = psm4_exdata$utility$lower,
                                             max = psm4_exdata$utility$upper),
                            dist = "unif")
psm_utility <- create_StateVals(utility_tbl, n = n_samples,
                                hesim_data = hesim_dat)

# Partitioned survival decision model
psm <- Psm$new(survival_models = psm_curves,
               utility_model = psm_utility,
               cost_models = list(medical = psm_costs_medical))
psm$sim_survival(t = seq(0, 5, 1/12))
autoplot(psm$survival_, labels = labs, ci = FALSE, ci_style = "ribbon")
psm$sim_stateprobs()
autoplot(psm$stateprobs_, labels = labs)
psm$sim_costs(dr = .03)
head(psm$costs_)
head(psm$sim_qalys(dr = .03)$qalys_)