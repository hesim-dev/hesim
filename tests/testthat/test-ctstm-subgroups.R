# Deterministic public-API model: unequal groups, distinct group outcomes,
# nonconsecutive IDs, two states with rewards, two strategies and PSA samples.
subgroup_fixture <- function(groups = c(90L, 10L, 10L, 10L)) {
  patients <- data.table(patient_id = c(11L, 35L, 60L, 92L), grp_id = groups)
  hd <- hesim_data(strategies = data.table(strategy_id = c(3L, 8L)),
                   patients = patients, states = data.table(state_id = 1:2))
  ids <- CJ(strategy_id = c(3L, 8L), patient_id = patients$patient_id,
            transition_id = 1:2)
  ids[patients, on = "patient_id", grp_id := i.grp_id]
  durations <- c(20, 2, 4, 6)[match(ids$patient_id, patients$patient_id)]
  durations[ids$transition_id == 2] <- 2
  durations <- durations * ifelse(ids$strategy_id == 3, 1, 1.2)
  inp <- input_mats(X = list(est = matrix(durations, ncol = 1)),
                    strategy_id = ids$strategy_id, patient_id = ids$patient_id,
                    grp_id = ids$grp_id, transition_id = ids$transition_id,
                    n_strategies = 2L, n_patients = 4L, n_transitions = 2L)
  par <- params_surv(coefs = list(est = matrix(c(1, 1.25), ncol = 1)),
                     dist = "fixed")
  tm <- matrix(NA_real_, 3, 3)
  tm[1, 2] <- 1
  tm[2, 3] <- 2
  trans <- IndivCtstmTrans$new(params = par, input_data = inp, trans_mat = tm)
  reward <- function(values) {
    create_StateVals(stateval_tbl(data.table(state_id = 1:2, est = values),
                                   dist = "fixed"), hesim_data = hd, n = 2)
  }
  model <- IndivCtstm$new(trans, reward(c(.8, .4)),
                          list(medical = reward(c(100, 300)),
                               other = reward(c(10, 0))))
  model$sim_disease(max_t = 10)

  # Analytic oracle, independent of the simulated trajectory and hesim's pv().
  truth <- CJ(sample = 1:2, strategy_id = c(3L, 8L),
               patient_id = patients$patient_id, state_id = 1:2, dr = c(0, .03))
  truth[patients, on = "patient_id", grp_id := i.grp_id]
  truth[, scale := ifelse(sample == 1, 1, 1.25) * ifelse(strategy_id == 3, 1, 1.2)]
  truth[, entry := c(20, 2, 4, 6)[match(patient_id, patients$patient_id)] * scale]
  truth[, start := ifelse(state_id == 1, 0, pmin(entry, 10))]
  truth[, stop := ifelse(state_id == 1, pmin(entry, 10), pmin(entry + 2 * scale, 10))]
  truth[, lys := ifelse(dr == 0, stop - start,
                         (exp(-dr * start) - exp(-dr * stop)) / dr)]
  truth[, qalys := lys * c(.8, .4)[state_id]]
  truth[, medical := lys * c(100, 300)[state_id]]
  truth[, other := lys * c(10, 0)[state_id]]
  list(model = model, truth = truth, patients = patients)
}

subgroup_expect <- function(actual, expected, keys, values) {
  expect_equal(nrow(actual), nrow(expected))
  a <- copy(actual)
  e <- copy(expected)
  setorderv(a, keys)
  setorderv(e, keys)
  for (col in c(keys, values)) expect_equal(a[[col]], e[[col]], tolerance = 1e-10)
}

test_that("CTSTM probabilities are conditional on subgroup", {
  f <- subgroup_fixture()
  m <- f$model
  m$sim_stateprobs(c(0, 3, 7, 10))
  p <- m$stateprobs_
  expect_equal(p[, sum(prob), by = .(sample, strategy_id, grp_id, t)]$V1,
               rep(1, 32))
  paths <- unique(f$truth[dr == 0, .(sample, strategy_id, patient_id, grp_id, entry, scale)])
  expected <- paths[, {
    times <- c(0, 3, 7, 10)
    state <- ifelse(times < entry, 1L, ifelse(times < entry + 2 * scale, 2L, 3L))
    data.table(t = rep(times, each = 3), state_id = rep(1:3, 4),
               prob = as.numeric(rep(1:3, 4) == rep(state, each = 3)))
  }, by = .(sample, strategy_id, patient_id, grp_id)]
  keys <- c("sample", "strategy_id", "grp_id", "state_id", "t")
  expected <- expected[, .(prob = mean(prob)), by = keys]
  subgroup_expect(p, expected, keys, "prob")
  # Labels must not depend on the first occurrence of a group in disprog.
  shuffled <- m$disprog_[order(-grp_id, -strategy_id)]
  subgroup_expect(m$trans_model$sim_stateprobs(t = c(0, 3, 7, 10), disprog = shuffled),
                  expected, keys, "prob")
})

test_that("CTSTM conditional means preserve population aggregation and membership", {
  for (groups in list(c(90L, 10L, 10L, 10L), rep(90L, 4))) {
    f <- subgroup_fixture(groups)
    m <- f$model
    for (by_patient in c(FALSE, TRUE)) {
      m$sim_qalys(dr = c(0, .03), by_patient = by_patient)
      m$sim_costs(dr = c(0, .03), by_patient = by_patient)
      keys <- c("sample", "strategy_id", "grp_id", "state_id", "dr")
      if (by_patient) keys <- c(keys, "patient_id")
      expected <- f$truth[, lapply(.SD, mean), by = keys,
                          .SDcols = c("qalys", "lys", "medical", "other")]
      subgroup_expect(m$qalys_, expected, keys, c("qalys", "lys"))
      for (cat in c("medical", "other")) {
        e <- copy(expected)
        e[, costs := get(cat)]
        subgroup_expect(m$costs_[category == cat], e, keys, "costs")
      }
      if (!by_patient) {
        expect_equal(m$qalys_$n_patients,
                     as.integer(table(groups)[as.character(m$qalys_$grp_id)]))
      } else {
        subgroup_expect(unique(m$qalys_[, .(patient_id, grp_id)]), f$patients,
                        "patient_id", "grp_id")
        expect_equal(m$qalys_[patient_id == 11 & state_id == 2, qalys], rep(0, 8))
      }
      for (by_grp in c(FALSE, TRUE)) {
        ce <- m$summarize(by_grp = by_grp)
        totals <- f$truth[, lapply(.SD, sum),
                          by = .(sample, strategy_id, patient_id, grp_id, dr),
                          .SDcols = c("qalys", "medical", "other")]
        keys <- c("sample", "strategy_id", "dr")
        if (by_grp) keys <- c(keys, "grp_id")
        expected <- totals[, lapply(.SD, mean), by = keys,
                           .SDcols = c("qalys", "medical", "other")]
        subgroup_expect(ce$qalys, expected, keys, "qalys")
        for (cat in c("medical", "other", "total")) {
          expected[, costs := if (cat == "total") medical + other else get(cat)]
          subgroup_expect(ce$costs[category == cat], expected, keys, "costs")
        }
        # Raw output is sufficient; no hidden model metadata is needed.
        direct <- hesim:::summarize_ce(copy(m$costs_), copy(m$qalys_), by_grp)
        setindexv(ce$costs, NULL)
        expect_equal(direct, ce)
      }
    }
  }
})

test_that("aggregate summaries retain legacy semantics and count each group once", {
  # Unequal numbers of state rows must not inflate the population denominator.
  q <- data.table(sample = 1L, strategy_id = 3L, dr = 0,
                  grp_id = c(90L, 10L, 10L), state_id = c(1L, 1L, 2L),
                  qalys = c(2, 1, 2))
  c <- copy(q)[, `:=`(costs = 100 * qalys, category = "medical", qalys = NULL)]
  legacy <- hesim:::summarize_ce(c, q)
  expect_equal(legacy$qalys$qalys, 5)
  expect_equal(legacy$costs$costs, c(500, 500))
  q[, n_patients := c(1L, 3L, 3L)]
  c[, n_patients := c(1L, 3L, 3L)]
  q[, qalys := qalys * 4 / n_patients]
  c[, costs := costs * 4 / n_patients]
  expect_equal(hesim:::summarize_ce(c, q), legacy)
  grouped <- hesim:::summarize_ce(c, q, by_grp = TRUE)
  expect_equal(grouped$qalys[order(grp_id), qalys], c(4, 8))
})
