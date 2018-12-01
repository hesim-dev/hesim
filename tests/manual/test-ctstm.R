context("Manual ctstm.R unit tests")
library("flexsurv")
library("mstate")
library("data.table")
rm(list = ls())

# Simulate flexsurv fit and compare state probabilities ------------------------
n_pats <- 10000
tmat <- rbind(c(NA, 1, 2), 
              c(NA, NA, 3), 
              c(NA, NA, NA))
transitions <- create_trans_dt(tmat)
transitions[, trans := factor(transition_id)]
hesim_dat <- hesim_data(strategies <- data.table(strategy_id = 1),
                        patients <- data.table(patient_id = 1:n_pats),
                        transitions = transitions)
hesim_edat <- expand(hesim_dat, by = c("strategies", "patients", "transitions"))

# Exponential model
fit <- flexsurvreg(Surv(years, status) ~ trans, data = bosms3, dist = "weibull")

## Simulate flexsurv
flexsurv_stateprobs <- pmatrix.simfs(fit, M = n_pats, t = 5, trans = tmat)

## Simulate hesim
transmod <- create_IndivCtstmTrans(fit, 
                                   data = hesim_edat, trans_mat = tmat,
                                   point_estimate = TRUE)
hesim_stateprobs <- transmod$sim_stateprobs(t = 5)

## Compare
stateprobs <- cbind(flexsurv_stateprobs[1, ],
                    hesim_stateprobs$prob)
colnames(stateprobs) <- c("flexsurv", "hesim")
print(stateprobs)
