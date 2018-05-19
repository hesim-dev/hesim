context("PartSurv.R unit tests")
library("flexsurv")
library("data.table")
library("pracma")

# Simulation
dt.strategies <- data.table(strategy_id = c(1, 2, 3))
dt.patients <- data.table(patient_id = seq(1, 3),
                          age = c(45, 50, 60),
                          female = c(0, 0, 1))
dt.states <- data.frame(state_id =  seq(1, 3),
                           state_name = paste0("state", seq(1, 3)))
hesim.dat <- hesim_data(strategies = dt.strategies,
                              patients = dt.patients,
                              states = dt.states)
N <- 5

# Partitioned survival curves  -------------------------------------------------
# Simulation data
curves.edata <- expand_hesim_data(hesim.dat, by = c("strategies", "patients"))

# Fit survival curves
surv.data <- part_surv4_simdata$survival
fits.exp <- fits.wei <- fits.splines <- vector(mode = "list", length = 3)
names(fits.exp) <- names(fits.wei) <- names(fits.splines) <- paste0("curves", seq(1, 3))
formulas <- list("Surv(endpoint1_time, endpoint1_status) ~ age",
                 "Surv(endpoint2_time, endpoint2_status) ~ age",
                 "Surv(endpoint3_time, endpoint3_status) ~ age")
for (i in 1:3){
  fits.exp[[i]] <- flexsurv::flexsurvreg(as.formula(formulas[[i]]),
                                         data = surv.data,
                                     dist = "exp")
  fits.wei[[i]] <- flexsurv::flexsurvreg(as.formula(formulas[[i]]), 
                                         data = surv.data,
                                          dist = "weibull")
  fits.splines[[i]] <- flexsurv::flexsurvspline(as.formula(formulas[[i]]), data = surv.data)
}
fits.exp <- partsurvfit(flexsurvreg_list(fits.exp), data = surv.data)
fits.wei <- partsurvfit(flexsurvreg_list(fits.wei), data = surv.data)
fits.splines <- partsurvfit(flexsurvreg_list(fits.splines), data = surv.data)

test_that("form_PartSurvCurves", {
  part.surv.curves <- form_PartSurvCurves(fits.wei, data = curves.edata, n = N,
                                          bootstrap = TRUE)
  expect_true(inherits(part.surv.curves, "PartSurvCurves"))
  expect_true(inherits(part.surv.curves$params, "params_surv_list"))
  expect_equal(as.numeric(part.surv.curves$data$X[[1]]$scale[, "age"]), 
              curves.edata$data$age)
  
  # errors
  expect_error(form_PartSurvCurves(3, data = curves.edata, n = N))
})

test_that("PartSurvCurves", {
  times <- c(1, 2, 3)
  
  # Sampling
  ## Weibull
  part.surv.curves <- form_PartSurvCurves(fits.wei, data = curves.edata, n = N)
  expect_true(inherits(part.surv.curves$survival(t = times), "data.table"))
  
  ## Splines
  part.surv.curves <- form_PartSurvCurves(fits.splines, data = curves.edata, n = N,
                                          bootstrap = FALSE)
  expect_equal(max(part.surv.curves$survival(t = times)$sample), N)
  
  # Comparison of summary of survival curves
  compare_surv_summary <- function(fits, data, fun_name = c("survival", "hazard",
                                                            "cumhazard", "rmst",
                                                            "quantile")){
    fun.name <- match.arg(fun_name)
    part.surv.curves <- form_PartSurvCurves(fits, data = data,
                                            point_estimate = TRUE,
                                            bootstrap = FALSE)
    hesim.out <- part.surv.curves[[fun.name]](t = times)
    fun.name2 <- if (fun.name == "cumhazard"){
      "cumhaz"
    } else{
      fun.name
    }
    flexsurv.out <- summary(fits$models[[1]], newdata = data.frame(age = data$data[1, age]),
                           t = times, type = fun.name2, tidy = TRUE, ci = FALSE)
    expect_equal(hesim.out[curve == 1, fun.name, with = FALSE][[1]], 
                 flexsurv.out[, "est"], tolerance = .001, scale = 1)
  }
  tmp.data <- curves.edata
  tmp.data$data <- tmp.data$data[1, ]
  
  compare_surv_summary(fits.wei, tmp.data, "survival")
  compare_surv_summary(fits.splines, tmp.data, "survival")
  
  compare_surv_summary(fits.wei, tmp.data, "hazard")
  compare_surv_summary(fits.splines, tmp.data, "hazard")
  
  compare_surv_summary(fits.wei, tmp.data, "cumhazard")
  compare_surv_summary(fits.splines, tmp.data, "cumhazard")
  
  compare_surv_summary(fits.wei, tmp.data, "rmst")
  compare_surv_summary(fits.splines, tmp.data, "rmst")
  
  # Quantiles
  part.surv.curves <- form_PartSurvCurves(fits.exp, data = curves.edata, n = N)
  X <- part.surv.curves$data$X$curves1$rate[1, , drop = FALSE]
  beta <- part.surv.curves$params$curves1$coefs$rate[1, , drop = FALSE]
  rate.hat <- X %*% t(beta)
  
  quantiles.out <- part.surv.curves$quantile(.5)
  expect_equal(qexp(.5, exp(rate.hat)), quantiles.out$quantile[1])
})

# Partitioned survival state values --------------------------------------------
fit.costs.medical <- stats::lm(costs ~ female + state_name, data = part_surv4_simdata$costs$medical)
edat <- expand_hesim_data(hesim.dat, by = c("strategies", "patients", "states"))
part.surv.costs.medical <- form_PartSurvStateVals(fit.costs.medical, data = edat, n = N)
part.surv.costs.medical2 <- form_PartSurvStateVals(fit.costs.medical, data = edat, n = N + 1)

test_that("PartSurvStateVals$predict", {
  expect_equal(c(part.surv.costs.medical$data$X %*% t(part.surv.costs.medical$params$coefs)),
              part.surv.costs.medical$predict()$value)
  
  expect_error(PartSurvStateVals$new(data = 3, params = 2)$predict())
  
  input.dat <- form_input_data(formula(~1), edat)
  expect_error(PartSurvStateVals$new(data = input.dat, params = 2)$predict())
})

# Partitioned survival model  --------------------------------------------------
set.seed(101)
times <- c(0, 2, 5, 8)
part.surv.curves <- form_PartSurvCurves(fits.wei, data = curves.edata, n = N)
part.surv.utility.data <- form_input_data(formula(~1), 
                                          expand_hesim_data(hesim.dat, 
                                                            by = c("strategies", "patients", "states")),
                                          id_vars = c("strategy_id", "patient_id", "state_id"))
part.surv.utility <- PartSurvStateVals$new(data = part.surv.utility.data,
                                           params = params_lm(coef = runif(N, .6, .8)))
part.surv <- PartSurv$new(survival_models = part.surv.curves,
                          utility_model = part.surv.utility,
                          cost_models = list(medical = part.surv.costs.medical))
expect_error(part.surv$sim_survival(t = c(2, 5)))
part.surv$sim_survival(t = times)

# State probabilities
test_that("PartSurv$stateprobs", {
  dt_by_grp <- function(x, by_var, value_var){
    df <- split(x, by = by_var)
    dt <- data.table(data.frame(lapply(df, function (x) x[[value_var]])))
  }
  
  surv.dt <- dt_by_grp(part.surv$survival_, by_var = "curve", value_var = "survival")
  surv.dt[, cross1 := ifelse(X1 > X2, 1, 0)]
  surv.dt[, cross2 := ifelse(X2 > X3, 1, 0)]
  n.crossings <- sum(surv.dt$cross1) + sum(surv.dt$cross2)
  if (n.crossings > 0){
    expect_warning(part.surv$sim_stateprobs()$stateprobs_)
  } else{
    part.surv$sim_stateprobs()$stateprobs_
  }
  
  state.probs.dt <- dt_by_grp(part.surv$stateprobs_, by_var = "state_id",
                              value_var = "prob")

  expect_equal(surv.dt$X1, state.probs.dt$X1)
  expect_equal(pmax(0, surv.dt$X2 - surv.dt$X1), state.probs.dt$X2)
  expect_equal(pmax(0, surv.dt$X3 - surv.dt$X2), state.probs.dt$X3)
  expect_equal(1 - surv.dt$X3, state.probs.dt$X4)
})

# Costs and QALYs
R_auc <- function(part_surv, type, type_num, dr = .03, 
                  state_id = 1, sample = 1, strategy_id = 1,
                  patient_id = 1){
  if (type == "costs_"){
    model <- part_surv$cost_models[[type_num]]
  } else{
    model <- part_surv$utility_model
  }
  dat <- model$data
  statevals <- dat$X %*% t(model$params$coefs)
  obs <- which(dat$state_id == state_id & dat$strategy_id == strategy_id &
                 dat$patient_id == patient_id)
  stateval <- statevals[obs, sample]
  
  env <- environment()
  stateprobs <- part.surv$stateprobs_[state_id == env$state_id &
                                      sample == env$sample &
                                      strategy_id == env$strategy_id &
                                      patient_id == env$patient_id]
  times <- stateprobs$t
  yvals <- exp(-dr * times) * stateval * stateprobs$prob
  return(pracma::trapz(x = times, y = yvals))
}

auc_compare <- function(part_surv, type = c("costs_", "qalys_"), 
                        type_num = NULL, dr,
                        state_id = 1, sample = 1, strategy_id = 1,
                  patient_id = 1){
  type <- match.arg(type)
  env <- environment()
  hesim.auc.dt <- part.surv[[type]][state_id == env$state_id &
                                  sample == env$sample &
                                  strategy_id == env$strategy_id &
                                  patient_id == env$patient_id &
                                  dr == env$dr]
  R.auc <- R_auc(part_surv, type = type, type_num = type_num, dr = dr, 
                 state_id = state_id, sample = sample,
                 strategy_id = strategy_id, patient_id = patient_id)
  expect_equal(hesim.auc.dt[[gsub("_", "", type)]], R.auc)
}

test_that("PartSurv$costs", {
  part.surv$sim_stateprobs()$stateprobs_
  part.surv$sim_costs(dr = c(0, .03))
  
  auc_compare(part.surv, type = "costs_", type_num = 1, dr = 0, strategy_id = 2)
  auc_compare(part.surv, type = "costs_", type_num = 1, dr = .03, strategy_id = 3)
  
  # Error messages
  part.surv2 <- PartSurv$new(survival_models = part.surv.curves,
                          utility_model = part.surv.utility,
                          cost_models = list(medical = part.surv.costs.medical2))
  part.surv2$sim_survival(t = times)
  expect_error(part.surv2$sim_costs(dr = c(0, .03)))
  part.surv2$sim_stateprobs()
  expect_error(part.surv2$sim_costs(dr = 0))
  
  ## Incorrect number of survival models
  partsurvfit2 <- partsurvfit(flexsurvreg_list(fits.wei$models[1:2]),
                              data = surv.data)
  part.surv.curves2 <- form_PartSurvCurves(partsurvfit2, 
                                           data = curves.edata, n = N,
                                          bootstrap = TRUE)
  part.surv2 <- PartSurv$new(survival_models = part.surv.curves2,
                             utility_model = part.surv.utility,
                             cost_models = list(medical = part.surv.costs.medical))
  part.surv2$sim_survival(t = times)
  part.surv2$sim_stateprobs()
  expect_error(part.surv2$sim_costs())
  
  ## Incorrect types
  part.surv2 <- PartSurv$new(survival_models = NULL)
  expect_error(part.surv2$sim_survival(t = times))
})

test_that("PartSurv$qalys", {
  part.surv$sim_stateprobs()$stateprobs_
  part.surv$sim_qalys(dr = c(0, .05))
  
  auc_compare(part.surv, type = "qalys_", dr = 0, strategy_id = 2)
  auc_compare(part.surv, type = "qalys_", dr = .05, patient_id = 2)
})