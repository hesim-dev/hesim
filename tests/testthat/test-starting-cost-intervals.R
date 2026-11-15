test_that("entry costs select the calendar interval independently for each entry", {
  hd <- hesim_data(strategies = data.frame(strategy_id = 1:2),
                   patients = data.frame(patient_id = 1L),
                   states = data.frame(state_id = 1:3))
  input <- expand(hd, by = c("strategies", "patients"))
  input[, intercept := 1]
  holding <- params_surv(list(est = matrix(2, dimnames = list(NULL, "intercept"))),
                         dist = "fixed")
  tm <- matrix(NA_real_, 4, 4)
  tm[cbind(1:3, 2:4)] <- 1:3
  trans <- create_IndivCtstmTrans(params_surv_list(holding, holding, holding),
                                  input, trans_mat = tm, clock = "reset")
  tariffs <- data.table::CJ(state_id = 1:3, time_start = c(0, 1, 10, 20, 30, 40))
  tariffs[, est := c(100, 1000, 10000, 20000, 30000, 40000)[match(time_start, c(0, 1, 10, 20, 30, 40))]]
  costs <- create_StateVals(stateval_tbl(tariffs, dist = "fixed"), n = 1,
                            hesim_data = hd, method = "starting")
  model <- IndivCtstm$new(trans_model = trans, cost_models = list(entry = costs))
  model$sim_disease(max_t = 7)
  for (dr in c(0, .03)) {
    model$sim_costs(dr = dr, by_patient = TRUE)
    expected <- c(100, 1000, 1000) * exp(-dr * c(0, 2, 4))
    expect_equal(model$costs_$costs, rep(expected, 2))
  }
  # Repeated entries must also remain in the terminal interval.
  tariffs <- tariffs[time_start <= 1]
  model$cost_models <- list(entry = create_StateVals(
    stateval_tbl(tariffs, dist = "fixed"), n = 1, hesim_data = hd, method = "starting"))
  model$sim_costs(dr = 0, by_patient = TRUE)
  expect_equal(model$costs_$costs, rep(c(100, 1000, 1000), 2))
})
