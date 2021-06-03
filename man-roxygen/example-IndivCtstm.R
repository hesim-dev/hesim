library("flexsurv")

# Treatment strategies, target population, and model structure
strategies <- data.frame(strategy_id = c(1, 2))
patients <- data.frame(patient_id = seq(1, 3),
                       age = c(45, 50, 60),
                       female = c(0, 0, 1))
states <- data.frame(state_id = c(1, 2))
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients,
                        states = states)

# Parameter estimation
## Multi-state model
tmat <- rbind(c(NA, 1, 2),
              c(3, NA, 4),
              c(NA, NA, NA))
fits <- vector(length = max(tmat, na.rm = TRUE), mode = "list")
surv_dat <- data.frame(mstate3_exdata$transitions)
for (i in 1:length(fits)){
  fits[[i]] <- flexsurvreg(Surv(years, status) ~ factor(strategy_id),
                           data = surv_dat,
                           subset = (trans == i),
                           dist = "weibull")
}
fits <- flexsurvreg_list(fits)

## Utility
utility_tbl <- stateval_tbl(data.frame(state_id = states$state_id,
                                       mean = mstate3_exdata$utility$mean,
                                       se = mstate3_exdata$utility$se),
                            dist = "beta")
## Costs
drugcost_tbl <- stateval_tbl(data.frame(strategy_id = strategies$strategy_id,
                                        est = mstate3_exdata$costs$drugs$costs),
                             dist = "fixed")
medcost_tbl <- stateval_tbl(data.frame(state_id = states$state_id,
                                       mean = mstate3_exdata$costs$medical$mean,
                                       se = mstate3_exdata$costs$medical$se),
                            dist = "gamma")

# Economic model
n_samples = 2

## Construct model
### Transitions
transmod_data <- expand(hesim_dat)
transmod <- create_IndivCtstmTrans(fits, input_data = transmod_data,
                                   trans_mat = tmat,
                                   n = n_samples)

### Utility
utilitymod <- create_StateVals(utility_tbl, n = n_samples, hesim_data = hesim_dat)

### Costs
drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples, hesim_data = hesim_dat)
medcostmod <- create_StateVals(medcost_tbl, n = n_samples, hesim_data = hesim_dat)
costmods <- list(drugs = drugcostmod,
                 medical = medcostmod)

### Combine
ictstm <- IndivCtstm$new(trans_model = transmod,
                         utility_model = utilitymod,
                         cost_models = costmods)


## Simulate outcomes
head(ictstm$sim_disease()$disprog_)
head(ictstm$sim_stateprobs(t = c(0, 5, 10))$stateprobs_[t == 5])
ictstm$sim_qalys(dr = .03)
ictstm$sim_costs(dr = .03)

### Summarize cost-effectiveness
ce <- ictstm$summarize()
head(ce)
format(summary(ce), pivot_from = "strategy")