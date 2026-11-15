test_that("linear-model predictions retain fitted factor encoding", {
  dat <- data.frame(group = factor(rep(c("A", "B"), each = 3)),
                    costs = c(90, 100, 110, 290, 300, 310))
  hd <- hesim_data(
    strategies = data.frame(strategy_id = 1L),
    patients = data.frame(patient_id = 1:2, group = factor(c("A", "B"))),
    states = data.frame(state_id = 1L)
  )
  input <- expand(hd, by = c("strategies", "patients", "states"))
  fits <- list(
    lm(costs ~ group, dat, contrasts = list(group = "contr.sum")),
    lm(costs ~ group, transform(dat, group = relevel(group, "B")))
  )
  for (fit in fits) {
    actual <- create_StateVals(fit, input, uncertainty = "none")$sim(0)$value
    expect_equal(actual, c(100, 300))
    expect_equal(actual, as.numeric(predict(fit, newdata = input)))
  }
})
