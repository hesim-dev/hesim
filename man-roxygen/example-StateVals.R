# Simple sick-sicker example where drug costs vary by treatment strategy
# and over time. Prior to time = 5, costs are $10,000 for treatment strategy 
# 1 and $5,000 for treatment strategy 2. After time = 5, costs are $2,000
# for both treatment strategies

## Setup the model
hesim_dat <- hesim_data(
  strategies = data.frame(strategy_id = c(1, 2)),
  patients <- data.frame(patient_id = 1:3),
  states = data.frame(state_id = c(1, 2), # Non-death states
                      state_name = c("sick", "sicker")) 
)

## Utility varies by health state
drugcost_tbl <- stateval_tbl(
  data.frame(
    strategy_id = c(1, 1, 2, 2),
    time_start = c(0, 5, 0, 5),
    est = c(10000, 2000, 5000, 2000)
  ),
  dist = "fixed"
)
drugcost_tbl

## Create drug cost model
drugcostmod <- create_StateVals(drugcost_tbl, n = 1, hesim_data = hesim_dat)
      
## Explore predictions from the drug cost model
drugcostmod$sim(t = c(2, 6), type = "predict")
 