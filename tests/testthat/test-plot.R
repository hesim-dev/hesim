context("plot.R unit tests")
library("data.table")

# Test autoplot method for stateprobs ------------------------------------------
# Create stateprobs object
trans <- rbind(
  c(5, 3, 2),
  c(1, 3, 6),
  c(0, 0, 1)
)
p <- array(NA, dim = c(3, 3, 4))
p[,, c(1, 3)] <- rdirichlet_mat(n = 2, trans)
p[,, c(2, 4)] <- apply_rr(p[,, c(1,3)], matrix(c(.5, .5), nrow = 2), index = list(c(1, 2)))
id <- data.table(sample = c(1, 1, 2, 2),
                 strategy_id = rep(c(1, 2), 2),
                 patient_id = 1)
tp <- tparams_transprobs(p, id)
mod <- CohortDtstmTrans$new(params = tp)
stprobs <- mod$sim_stateprobs(n_cycles = 20)

labs <- list("strategy_id" = c("s1" = 1, 
                               "s2" = 2),
             "state_id" = c("state1" = 1, 
                            "state2" = 2,
                            "state3" = 3))

test_that("autoplot.stateprobs() returns ggplot", {
  expect_true(inherits(autoplot(stprobs), "ggplot"))
})

test_that("autoplot.stateprobs() correctly passes labels", {
  p <- autoplot(stprobs, labels = labs)
  expect_equal(levels(p$data$strategy_id), names(labs$strategy_id))
  expect_equal(levels(p$data$state_id), names(labs$state_id))
})

test_that("autoplot.stateprobs() allows confidence intervals", {
  # Confidence intervals as lines
  p <- autoplot(stprobs, labels = labs, ci = TRUE)
  expect_equal(p$labels$fill, "strategy_id")
  
  # Confidence intervals as ribbon
  p <- autoplot(stprobs, labels = labs, ci = TRUE, ci_style = "ribbon")
  expect_equal(p$labels$fill, "strategy_id")

})

test_that("autoplot.stateprobs() works with patient weights", {
  stprobs2 <- copy(stprobs)[, c("patient_id", "grp_id") := 2L]
  stprobs2 <- rbind(stprobs, stprobs2)
  stprobs2[, patient_wt := ifelse(patient_id == 1, .3, .7)]
  setattr(stprobs2, "class", c("stateprobs", "data.table", "data.frame"))
  expect_true(inherits(autoplot(stprobs2), "ggplot"))
})

