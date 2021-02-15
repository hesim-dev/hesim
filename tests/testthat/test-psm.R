context("psm.R unit tests")
library("flexsurv")
library("data.table")
library("pracma")
rm(list = ls())

# hesim data
strategies_dt <- data.table(strategy_id = seq(2, 4)) # testing for cases when doesn't start at 1
patients_dt <- data.table(patient_id = seq(1, 3),
                          patient_wt = 1,
                          age = c(45, 50, 60),
                          female = c(0, 0, 1))
states_dt <- data.frame(state_id =  seq(1, 3),
                        state_name = paste0("state", seq(1, 3)))
hesim_dat <- hesim_data(strategies = strategies_dt,
                        patients = patients_dt,
                        states = states_dt)
N <- 5

# Partitioned survival curves  -------------------------------------------------
# Simulation data
surv_input_data <- expand(hesim_dat, by = c("strategies", "patients"))

# Fit survival curves
surv_est_data <- psm4_exdata$survival
fits_exp <- fits_wei <- fits_weinma <- fits_spline <- fits_ggamma <- vector(mode = "list", length = 3)
names(fits_exp) <- names(fits_wei) <- names(fits_spline) <- paste0("curves", seq(1, 3))
formulas <- list("Surv(endpoint1_time, endpoint1_status) ~ age",
                 "Surv(endpoint2_time, endpoint2_status) ~ age",
                 "Surv(endpoint3_time, endpoint3_status) ~ age")
for (i in 1:3){
  fits_exp[[i]] <- flexsurv::flexsurvreg(as.formula(formulas[[i]]),
                                         data = surv_est_data,
                                        dist = "exp")
  fits_wei[[i]] <- flexsurv::flexsurvreg(as.formula(formulas[[i]]), 
                                         data = surv_est_data,
                                          dist = "weibull")
  fits_weinma[[i]] <- suppressWarnings(flexsurv::flexsurvreg(as.formula(formulas[[i]]), 
                                         data = surv_est_data,
                                         dist = hesim_survdists$weibullNMA,
                                         inits = fits_wei[[i]]$res.t[, "est"]))
  fits_spline[[i]] <- flexsurv::flexsurvspline(as.formula(formulas[[i]]), data = surv_est_data)
  fits_ggamma[[i]] <- flexsurv::flexsurvreg(as.formula(formulas[[i]]),
                                         data = surv_est_data,
                                        dist = "gengamma")
}
fits_exp <- flexsurvreg_list(fits_exp)
fits_wei <- flexsurvreg_list(fits_wei)
fits_weinma <- flexsurvreg_list(fits_weinma)
fits_spline <- flexsurvreg_list(fits_spline)
fits_ggamma <- flexsurvreg_list(fits_ggamma)

test_that("create_PsmCurves() produces expected output", {
  psm_curves <- create_PsmCurves(fits_wei, input_data = surv_input_data, n = N,
                                 uncertainty = "bootstrap", est_data = surv_est_data)
  expect_true(inherits(psm_curves, "PsmCurves"))
  expect_true(inherits(psm_curves$params, "params_surv_list"))
  expect_equal(as.numeric(psm_curves$input_data$X[[1]]$scale[, "age"]), 
              surv_input_data$age)
  expect_equal(colnames(psm_curves$params$curves1$coefs$shape)[1],
               "(Intercept)")
})

test_that("create_PsmCurves() returns errors when expected", {
  expect_error(
    create_PsmCurves(fits_wei, input_data = surv_input_data, n = N,
                    uncertainty = "bootstrap"),
    "If uncertainty == 'bootstrap', then est_data cannot be NULL"
  )
})

test_that("PsmCurves are correct", {
  times <- c(1, 2, 3)
  
  # Sampling
  ## Weibull
  psm_curves <- create_PsmCurves(fits_wei, input_data = surv_input_data, n = N,
                                 uncertainty = "bootstrap",
                                 est_data = surv_est_data)
  expect_true(inherits(psm_curves$survival(t = times), "data.table"))
  
  ## Splines
  psm_curves <- create_PsmCurves(fits_spline, input_data = surv_input_data, n = N,
                                uncertainty = "normal")
  expect_equal(max(psm_curves$survival(t = times)$sample), N)
  
  # Comparison of summary of survival curves
  compare_surv_summary <- function(fits, data, fun_name = c("survival", "hazard",
                                                            "cumhazard", "rmst",
                                                            "quantile")){
    fun_name <- match.arg(fun_name)
    psm_curves <- create_PsmCurves(fits, input_data = data,
                                   uncertainty = "none")
    
    hesim_out <- psm_curves[[fun_name]](t = times)
    fun_name2 <- if (fun_name == "cumhazard"){
      "cumhaz"
    } else{
      fun_name
    }
    flexsurv_out <- summary(fits[[1]], newdata = data.frame(age = data[1, age]),
                           t = times, type = fun_name2, tidy = TRUE, ci = FALSE)
    expect_equal(hesim_out[curve == 1, fun_name, with = FALSE][[1]], 
                 flexsurv_out[, "est"], tolerance = .001, scale = 1)
  }
  tmp_data = surv_input_data
  tmp_data <- tmp_data[1, ]
  
  compare_surv_summary(fits_wei, tmp_data, "survival")
  compare_surv_summary(fits_spline, tmp_data, "survival")
  compare_surv_summary(fits_weinma, tmp_data, "survival")
  compare_surv_summary(fits_ggamma, tmp_data, "survival")
  
  compare_surv_summary(fits_wei, tmp_data, "hazard")
  compare_surv_summary(fits_spline, tmp_data, "hazard")
  compare_surv_summary(fits_weinma, tmp_data, "hazard")
  compare_surv_summary(fits_ggamma, tmp_data, "hazard")
  
  
  compare_surv_summary(fits_wei, tmp_data, "cumhazard")
  compare_surv_summary(fits_spline, tmp_data, "cumhazard")
  compare_surv_summary(fits_weinma, tmp_data, "cumhazard")
  compare_surv_summary(fits_ggamma, tmp_data, "cumhazard")
  
  compare_surv_summary(fits_wei, tmp_data, "rmst")
  compare_surv_summary(fits_spline, tmp_data, "rmst")
  compare_surv_summary(fits_weinma, tmp_data, "rmst")
  
  # Quantiles
  psm_curves <- create_PsmCurves(fits_exp, input_data = surv_input_data, n = N,
                                 uncertainty = "bootstrap", est_data = surv_est_data)
  X <- psm_curves$input_data$X$curves1$rate[1, , drop = FALSE]
  beta <- psm_curves$params$curves1$coefs$rate[1, , drop = FALSE]
  rate_hat <- X %*% t(beta)
  
  quantiles_out <- psm_curves$quantile(.5)
  expect_equal(qexp(.5, exp(rate_hat)), quantiles_out$quantile[1])
})

# Partitioned survival model  --------------------------------------------------
set.seed(101)

# Construct PSM
times <- c(0, 2, 5, 8)

## Survival models
psm_curves <- create_PsmCurves(fits_wei, input_data = surv_input_data, n = N,
                               uncertainty = "normal")

## Utility model
psm_X <- create_input_mats(formula_list(mu = formula(~1)), 
                                     expand(hesim_dat, 
                                     by = c("strategies", "patients", "states")),
                                     id_vars = c("strategy_id", "patient_id", "state_id"))
psm_utility <- StateVals$new(input_data = psm_X,
                             params = params_lm(coef = runif(N, .6, .8)))


## Cost model(s)
fit_costs_medical <- stats::lm(costs ~ female + state_name, 
                               data = psm4_exdata$costs$medical)
cost_input_data <- expand(hesim_dat, by = c("strategies", "patients", "states"))
psm_costs_medical <- create_StateVals(fit_costs_medical, 
                                      input_data = cost_input_data, 
                                      n = N)
psm_costs_medical2 <- create_StateVals(fit_costs_medical, 
                                       input_data = cost_input_data, 
                                       n = N + 1)

## Combine
psm <- Psm$new(survival_models = psm_curves,
               utility_model = psm_utility,
               cost_models = list(medical = psm_costs_medical))
psm$sim_survival(t = times)

# Run tests
## $sim_survival()
test_that("$sim_survival() returns an error if the first element of t is not 0", {
  expect_error(psm$sim_survival(t = c(2, 5)),
               "The first element of 't' must be 0") 
})

## $sim_stateprobs()
dt_by_grp <- function(x, by_var, value_var){
  df <- split(x, by = by_var)
  dt <- data.table(data.frame(lapply(df, function (x) x[[value_var]])))
}

surv_dt <- dt_by_grp(psm$survival_, by_var = "curve", value_var = "survival")
surv_dt[, cross1 := ifelse(X1 > X2, 1, 0)]
surv_dt[, cross2 := ifelse(X2 > X3, 1, 0)]
n_crossings <- sum(surv_dt$cross1) + sum(surv_dt$cross2)

test_that("$sim_stateprobs() produces warning if curves cross", {
  if (n_crossings > 0){
    expect_warning(psm$sim_stateprobs()$stateprobs_)
  } else{ # No warning expected
    expect_warning(psm$sim_stateprobs()$stateprobs_, NA)
  }
})

test_that("$sim_stateprobs() should return patient_wt if not NULL in hesim_data", {
  psm$sim_stateprobs()
  expect_true(!is.null(psm$stateprobs_$patient_wt))
})

test_that("sim_stateprobs produces expected results", {
  psm$sim_stateprobs()
  stateprobs_dt <- dt_by_grp(psm$stateprobs_, by_var = "state_id",
                              value_var = "prob")

  expect_equal(surv_dt$X1, stateprobs_dt$X1)
  expect_equal(pmax(0, surv_dt$X2 - surv_dt$X1), stateprobs_dt$X2)
  expect_equal(pmax(0, surv_dt$X3 - surv_dt$X2), stateprobs_dt$X3)
  expect_equal(1 - surv_dt$X3, stateprobs_dt$X4)
})

## sim_costs() and sim_qalys()
test_that("sim_costs() and sim_qalys() both return a data.table", {
  psm$sim_stateprobs()
  psm$sim_costs(dr = c(0, .03))
  expect_true(inherits(psm$costs_, "data.table"))
  psm$sim_qalys(dr = c(0, .05))
  expect_true(inherits(psm$qalys_, "data.table"))
})

## Psm from a parameter object
test_that("A Psm object can be constructed and simulated from a parameter object", {
  # PsmCurves
  params_wei <- create_params(fits_wei, n = 5)
  tmp_input_data <- surv_input_data
  tmp_input_data[["(Intercept)"]] <- 1
  psm_curves <- create_PsmCurves(params_wei, input_data = tmp_input_data)
  expect_true(inherits(psm_curves$hazard(t = c(1, 2, 3)),
                      "data.table"))
  
  # Psm
  psm <- Psm$new(survival_models = psm_curves)
  psm$sim_survival(t = c(0, 1, 2, 3))
  expect_true(inherits(psm$survival_,
                      "data.table")) 
  psm$sim_stateprobs()
  expect_true(inherits(psm$stateprobs_,
                      "data.table")) 
})
