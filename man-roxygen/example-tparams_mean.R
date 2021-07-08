# Setup model
hesim_dat <- hesim_data(
  strategies = data.frame(strategy_id = c(1, 2)),
  patients = data.frame(patient_id = c(1, 2)),
  states = data.frame(
    state_id = c(1, 2, 3),
    state_name = c("state1", "state2", "state3")
  )
)

# Cost model
cost_tbl <- stateval_tbl(
  data.frame(strategy_id = hesim_dat$strategies$strategy_id,
             mean = c(5000, 3000),
             se = c(200, 100)
            ),
  dist = "gamma"
)
costmod <- create_StateVals(cost_tbl, n = 2, hesim_data = hesim_dat)

# The 'params' field of the `StateVals` class is a tparams_mean object
class(costmod$params)
costmod$params
summary(costmod$params)
