test_that("custom state values follow identifiers, not row position", {
  x <- data.table::CJ(sample = 1:2, strategy_id = 1:2, patient_id = 1:2,
                      state_id = 1:2, time_start = c(0, 5))
  x[, value := sample + 10 * strategy_id + 100 * patient_id +
                 1000 * state_id + 10000 * time_start]
  # Reverse only the second draw, so within-draw orders disagree.
  shuffled <- rbind(x[sample == 1], x[sample == 2][.N:1])
  out <- create_StateVals(stateval_tbl(shuffled, dist = "custom"), n = 2)$sim(c(0, 6))
  ids <- c("sample", "strategy_id", "patient_id", "state_id", "time")
  expected <- data.table::copy(x)
  data.table::setnames(expected, "time_start", "time")
  expected[time == 5, time := 6]
  data.table::setorderv(out, ids)
  data.table::setorderv(expected, ids)
  expect_equal(out$value, expected$value)

  # Group-level tables must align before expanding to individual patients.
  hd <- hesim_data(strategies = data.frame(strategy_id = 1L),
                   patients = data.frame(patient_id = 1:2, risk = c("low", "high")),
                   states = data.frame(state_id = 1L))
  groups <- data.frame(sample = c(1L, 1L, 2L, 2L),
                       risk = c("low", "high", "high", "low"),
                       state_id = 1L, value = c(.9, .5, .4, .8))
  out <- create_StateVals(stateval_tbl(groups, dist = "custom", grp_var = "risk"),
                          hesim_data = hd, n = 2)$sim(0)
  expect_equal(out$value, c(.9, .5, .8, .4))
})
