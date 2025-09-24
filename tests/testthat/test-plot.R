context("plot.R unit tests")
library("data.table")

get_labs <- function(x) x$labels
if ("get_labs" %in% getNamespaceExports("ggplot2")) {
  get_labs <- ggplot2::get_labs
}

# Test autoplot method for survival --------------------------------------------
# Create mock survival object
hesim_dat <- hesim_data(
  strategies = data.frame(strategy_id = 1:2L, 
                          strategy2 = c(0, 1)),
  patients = data.frame(patient_id = 1L, cons = 1)
)
params1 <- params_surv(
  coefs = list(rate = as.matrix(data.frame(
    cons = c(log(.2), log(.25)),
    strategy2 = c(log(.8), log(.9))
    ))),
  dist = "exp"
)
params2 <- params_surv(
  coefs = list(rate = as.matrix(data.frame(
    cons = c(log(.1), log(.15)),
    strategy2 = c(log(.8), log(.9))
    ))),
  dist = "exp"
)
params <- params_surv_list(params1, params2)
psm_curves <- create_PsmCurves(params, input_data = expand(hesim_dat))
surv <- psm_curves$survival(t = 0:20)

labs <- list(
  "strategy_id" = c("s1" = 1, 
                    "s2" = 2),
  "curve" = c("PFS" = 1, 
              "OS" = 2),
  "state_id" = c("state1" = 1, 
                 "state2" = 2,
                  "state3" = 3)
)

test_that("autoplot.survival() returns ggplot", {
  expect_true(inherits(autoplot(surv), "ggplot"))
})

test_that("autoplot.survival() correctly passes labels", {
  p <- autoplot(surv, labels = labs)
  expect_equal(levels(p$data$strategy_id), names(labs$strategy_id))
  expect_equal(levels(p$data$curve), names(labs$curve))
})

test_that("autoplot.survival() allows confidence intervals", {
  # Confidence intervals as lines
  p <- autoplot(surv, labels = labs, ci = TRUE)
  
  expect_equal(get_labs(p)$colour, "Curve")
  # Confidence intervals as ribbon
  p <- autoplot(surv, labels = labs, ci = TRUE, ci_style = "ribbon")
  expect_equal(get_labs(p)$colour, "Curve")
})

test_that("autoplot.survival() works with patient weights", {
  surv[, patient_wt := ifelse(patient_id == 1, .3, .7)]
  expect_true(inherits(autoplot(surv), "ggplot"))
})

# Test autoplot method for stateprobs ------------------------------------------
# Create mock stateprobs object
trans <- rbind(
  c(5, 3, 2),
  c(1, 3, 6),
  c(0, 0, 1)
)
p <- array(NA, dim = c(3, 3, 4))
p[,, c(1, 3)] <- rdirichlet_mat(n = 2, trans)
p[,, c(2, 4)] <- apply_rr(p[,, c(1,3)], matrix(c(.5, .5), nrow = 2), 
                          index = list(c(1, 2)))
id <- data.table(sample = c(1, 1, 2, 2),
                 strategy_id = rep(c(1, 2), 2),
                 patient_id = 1)
tp <- tparams_transprobs(p, id)
mod <- CohortDtstmTrans$new(params = tp)
stprobs <- mod$sim_stateprobs(n_cycles = 20)

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
  
  expect_equal(get_labs(p)$colour, "Strategy")
  # Confidence intervals as ribbon
  p <- autoplot(stprobs, labels = labs, ci = TRUE, ci_style = "ribbon")
  expect_equal(get_labs(p)$colour, "Strategy")

})

test_that("autoplot.stateprobs() works with patient weights", {
  stprobs2 <- copy(stprobs)[, c("patient_id", "grp_id") := 2L]
  stprobs2 <- rbind(stprobs, stprobs2)
  stprobs2[, patient_wt := ifelse(patient_id == 1, .3, .7)]
  setattr(stprobs2, "class", c("stateprobs", "data.table", "data.frame"))
  expect_true(inherits(autoplot(stprobs2), "ggplot"))
})
