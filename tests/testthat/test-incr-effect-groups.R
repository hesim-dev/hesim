test_that("incremental effects pair the same subgroup and PSA draw", {
  x <- data.table::CJ(group = c("high", "low"), arm = 1:3, draw = 1:3)
  x[, costs := ifelse(group == "high", 10000, 1000) + 100 * draw + 50 * arm]
  x[, qalys := ifelse(group == "high", 5, 10) + 0.1 * draw + 0.5 * arm]
  # The comparator is arm 1: arm 2 adds 50/0.5; arm 3 adds 100/1.
  for (arms in list(1:3, 1:2)) {
    input <- x[arm %in% arms]
    input <- input[rev(seq_len(nrow(input)))]
    out <- incr_effect(input, comparator = 1, sample = "draw", strategy = "arm",
                       grp = "group", outcomes = c("costs", "qalys"))
    expect_equal(nrow(out), 6L * (length(arms) - 1L))
    expect_equal(out$icosts, 50 * (out$arm - 1))
    expect_equal(out$iqalys, 0.5 * (out$arm - 1))
  }
})
