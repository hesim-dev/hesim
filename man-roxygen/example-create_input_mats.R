library("flexsurv")

strategies <- data.frame(strategy_id = c(1, 2))
patients <- data.frame(patient_id = seq(1, 3),
                       age = c(45, 47, 60),
                       female = c(1, 0, 0),
                       group = factor(c("Good", "Medium", "Poor")))
states <- data.frame(state_id =  seq(1, 3),
                     state_name = factor(paste0("state", seq(1, 3))))
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients,
                        states = states)

# Class "lm"
input_data <- expand(hesim_dat, by = c("strategies", "patients", "states"))
medcost_fit <- lm(costs ~ female + state_name, psm4_exdata$costs$medical)
input_mats <- create_input_mats(fit_lm, input_data)
input_mats

# Class "flexsurvreg"
input_data <- expand(hesim_dat, by = c("strategies", "patients"))
fit_wei <- flexsurvreg(formula = Surv(futime, fustat) ~ 1,
                       data = ovarian, dist = "weibull")
input_mats <- create_input_mats(fit_wei, input_data)
input_mats