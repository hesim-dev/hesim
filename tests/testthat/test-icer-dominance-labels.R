test_that("dominance follows strategy and subgroup label order", {
  x <- data.table::CJ(sample = 1:2, strategy = 1:3, group = 1:2)
  x[, qalys := 5 + c(0, 1, -1)[strategy] * ifelse(group == 1, 1, -1)]
  x[, costs := 20000 + c(0, -10000, 10000)[strategy] * ifelse(group == 1, 1, -1)]
  pw <- cea_pw(x, comparator = 1, sample = "sample", strategy = "strategy",
               grp = "group", e = "qalys", c = "costs", k = 50000)
  for (arms in list(c(1, 2, 3), c(3, 2, 1))) {
    for (groups in list(c(1, 2), c(2, 1))) {
      labels <- list(strategy = setNames(arms, paste0("arm", arms)),
                     group = setNames(groups, paste0("group", groups)))
      out <- format(icer(pw, labels = labels), pivot_from = NULL,
                     pretty_names = FALSE, drop_grp = FALSE)
      out <- out[outcome == "ICER"]
      expect_identical(out[strategy == "arm2" & grp == "group1", value], "Dominates")
      expect_identical(out[strategy == "arm3" & grp == "group1", value], "Dominated")
      expect_identical(out[strategy == "arm2" & grp == "group2", value], "Dominated")
      expect_identical(out[strategy == "arm3" & grp == "group2", value], "Dominates")
    }
  }
})
