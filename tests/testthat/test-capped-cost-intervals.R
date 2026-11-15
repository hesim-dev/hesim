test_that("capped courses integrate only their overlap with tariff intervals", {
  hd <- hesim_data(strategies = data.frame(strategy_id = 1:2),
                   patients = data.frame(patient_id = 1L, intercept = 1),
                   states = data.frame(state_id = 1:2))
  holding <- function(t) params_surv(list(est = data.frame(intercept = t)), dist = "fixed")
  tm <- rbind(c(NA, 1, NA), c(NA, NA, 2), c(NA, NA, NA))
  trans <- create_IndivCtstmTrans(params_surv_list(holding(3), holding(2)),
                                  expand(hd), trans_mat = tm, clock = "reset")
  tariffs <- data.frame(state_id = rep(1:2, each = 2),
                        time_start = rep(c(0, 2), 2), est = rep(c(10000, 2000), 2))
  costs <- create_StateVals(stateval_tbl(tariffs, dist = "fixed"), n = 1, hesim_data = hd)
  model <- IndivCtstm$new(trans_model = trans, cost_models = list(course = costs))
  model$sim_disease(max_t = 10)
  # Episodes are [0,3] and [3,5]. Integrate intersections independently.
  for (reset in c(FALSE, TRUE)) for (cap in c(0, 1, 2.5, Inf)) for (dr in c(0, .03)) {
    costs$time_reset <- reset
    model$sim_costs(dr = dr, max_t = cap, by_patient = TRUE)
    expected <- vapply(1:2, function(h) {
      start <- c(0, 3)[h]
      stop <- min(c(3, 5)[h], start + cap)
      shift <- if (reset) start else 0
      a <- pmax(start, c(0, 2) + shift)
      b <- pmin(stop, c(2, Inf) + shift)
      keep <- b > a
      if (dr == 0) sum(c(10000, 2000)[keep] * (b[keep] - a[keep])) else
        sum(c(10000, 2000)[keep] * (exp(-dr * a[keep]) - exp(-dr * b[keep])) / dr)
    }, numeric(1))
    expect_equal(model$costs_$costs, rep(expected, 2))
    expect_true(all(model$costs_$costs >= 0))
  }
})
