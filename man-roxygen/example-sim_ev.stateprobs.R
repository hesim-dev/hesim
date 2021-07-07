# We need (i) a state probability object and (ii) a model for state values
## We should start by setting up our decision problem
hesim_dat <-  hesim_data(strategies = data.frame(strategy_id = 1:2),
                         patients = data.frame(patient_id = 1:3),
                         states = data.frame(state_id = 1))
input_data <- expand(hesim_dat, by = c("strategies", "patients"))

## (i) Simulate a state probability object
tpmat_id <- tpmatrix_id(input_data, n_samples = 2) 
p_12 <- ifelse(tpmat_id$strategy_id == 1, .15, .1)
tpmat <- tpmatrix(
  C, p_12,
  0, 1
)
transmod <- CohortDtstmTrans$new(params = tparams_transprobs(tpmat, tpmat_id))
stprobs <- transmod$sim_stateprobs(n_cycles = 3)

## Construct model for state values
outcome_tbl <- stateval_tbl(
  data.frame(
    state_id = 1,
    est = 5000
  ),
  dist = "fixed"
)
outmod <- create_StateVals(outcome_tbl, n = 2, hesim_data = hesim_dat)

# We can then simulate expected values
## The generic expected values function
sim_ev(stprobs, models = outmod)

## We can also pass a list of models
sim_ev(stprobs, models = list(`Outcome 1` = outmod))

## Suppose the outcome were a cost category. Then we might
## prefer the following:
sim_costs(stprobs, models = list(drug = outmod))

## Length of stay is computed if there is no state value model
sim_ev(stprobs)
