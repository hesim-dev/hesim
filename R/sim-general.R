# Disease progression object ---------------------------------------------------
#' Disease progression object
#'
#' An object of class `disprog` returned from methods 
#' `$sim_disease()` in model classes. It contains simulated trajectories 
#' through a multi-state model.  
#' 
#' @section Components:
#' A `disprog` object inherits from `data.table` and contains
#' the following columns:
#' 
#' \describe{
#' \item{sample}{A random sample from the PSA.}
#' \item{strategy_id}{The treatment strategy ID.}
#' \item{patient_id}{The patient ID.}
#' \item{from}{The health state ID transitioned from.}
#' \item{to}{The health state ID transitioned to.}
#' \item{final}{An indicator equal to 1 if a patient is in their final health
#' state during the simulation and 0 otherwise.}
#' \item{time_start}{The time at the start of the interval.}
#' \item{time_stop}{The time at the end of the interval.}
#' }
#'
#' @seealso [`IndivCtstm`], [`IndivCtstmTrans`]
#' @name disprog
NULL

# Survival ---------------------------------------------------------------------
#' Survival object
#'
#' An object of class `survival` stores survival probabilities. It is typically
#' returned by `Psm$sim_survival()` or `PsmCurves$survival()`; however, it can also
#' be constructed "manually" from existing data using the `survival()` 
#' function as described below. The latter option is useful if survival modeling
#' has been performed by an `R` package other than those that integrate with `hesim` (
#' currently `flexsurv`). In this case a simulation model can still be developed
#' by using [`sim_stateprobs.survival()`] to compute simulated state probabilities and
#' then simulating quality-adjusted life-years and costs in a typical fashion.
#' 
#' @param data A tabular object that can be coerced to a `data.table` with
#' [`as.data.table()`].
#' @param sample The name of the column corresponding to `sample`.
#' @param strategy_id The name of the column corresponding to `strategy_id`.
#' @param patient_id The name of the column corresponding to `patient_id`.
#' @param grp_id The name of the column corresponding to `grp_id`.
#' @param curve The name of the column corresponding to `curve`.
#' @param t The name of the column corresponding to `t`.
#' @param survival The name of the column corresponding to `survival`.
#' 
#' @return
#' An object of class `survival` that inherits from `data.table` and contains
#' the following columns:
#' 
#' \describe{
#'   \item{sample}{A random sample from the PSA.}
#'   \item{strategy_id}{The treatment strategy ID.}
#'   \item{patient_id}{The patient ID.}
#'   \item{grp_id}{The subgroup ID.}
#'   \item{curve}{One of the `N`-1 survival curves in an N-state partitioned
#'   survival model. Each curve corresponds to unique endpoint.}
#'   \item{t}{The time at which a survival probability is computed.}
#'   \item{survival}{The probability of surviving to time `t`.}
#' }
#' The object also contains a `size` attribute that contains the elements
#' `n_samples`, `n_strategies`, `n_patients`, `n_states`, and `n_times` denoting
#'  the number of samples, treatment strategies, patients, health states, and times. 
#'
#' @seealso `survival` objects are returned by methods in the [`Psm`] and [`PsmCurves`] 
#' classes. An example in which a `survival` object is constructed "manually"
#' (presumably from a preexisting survival model fit using software other than `flexsurv`) 
#' is provided in the documentation to [`sim_stateprobs.survival()`].
#' @export
survival <- function(data, sample = "sample", strategy_id = "strategy_id",
                     patient_id = "patient_id", grp_id = "grp_id",
                     curve = "curve", t = "t", survival = "survival") {
  x <- as.data.table(data)
  
  # Get right columns
  user_cols <- c(sample, strategy_id, patient_id, grp_id, curve, t, survival)
  cols <- c("sample", "strategy_id", "patient_id", "grp_id", "curve", "t", "survival")
  x <- x[, user_cols, with = FALSE]
  setnames(x, user_cols, cols)
  
  # Make sure sorted correctly
  setorderv(x, c("sample", "strategy_id", "patient_id", "grp_id", "curve", "t"))
  
  # Some checks
  ## Number of total observations is correct
  get_n <- function(v) length(unique(x[[v]]))
  id_cols <- c("sample", "strategy_id", "patient_id", "grp_id", "curve", "t")
  id_n <- sapply(id_cols, get_n)
  if(prod(id_n) != nrow(x)) {
    stop(paste0("The number of rows in 'data' must be equal to the product of the ", 
                "number of unique values of the 'sample', 'strategy_id', 'patient_id' ",
                "'grp_id', 'curve', and 't' columns."))
  }
  
  # Return
  setattr(x, "class", c("survival", "data.table", "data.frame"))
  setattr(x, "size", 
          c(n_samples = id_n[["sample"]],
            n_strategies = id_n[["strategy_id"]],
            n_patients = id_n[["patient_id"]],
            n_states = id_n[["curve"]] + 1,
            n_times = id_n[["t"]]))
  return(x)
}

# State probabilities ----------------------------------------------------------
#' State probability object
#'
#' An object of class `stateprobs` returned by [`sim_stateprobs()`] or from 
#' `$sim_stateprobs()` methods in model classes. 
#' 
#' @section Components:
#' A `stateprobs` object inherits from `data.table` and contains
#' the following columns:
#' 
#' \describe{
#'   \item{sample}{A random sample from the PSA.}
#'   \item{strategy_id}{The treatment strategy ID.}
#'   \item{patient_id}{The patient ID.}
#'   \item{grp_id}{The subgroup ID.}
#'   \item{state_id}{The health state ID.}
#'   \item{t}{The time at which a state probability is computed.}
#'   \item{prob}{The probability of being in a given health state.}
#' }
#' 
#' When simulating individual-level models, the `patient_id` column is
#' not included as state probabilities are computed by averaging across patients.
#' 
#' In cohort models, the object also contains a `size` attribute that contains 
#' the elements `n_samples`, `n_strategies`, `n_patients`, `n_states`, and
#'  `n_times` denoting the number of samples, treatment strategies, patients, 
#'  health states, and times. 
#'
#' @name stateprobs
NULL

#' Simulated state probabilities
#'
#' A generic function to simulate state probabilities and create an object of 
#' class [`stateprobs`].
#' 
#' @param x An object of the appropriate class.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' A [`stateprobs`] object.
#' @export
#' @seealso [`sim_stateprobs.survival`]
sim_stateprobs <- function(x, ...) {
  UseMethod("sim_stateprobs")
}

#' Simulate state probabilities from survival curves
#'
#' Simulate health state probabilities from a [`survival`] object using partitioned
#' survival analysis.
#' 
#' @param x An object of class [`survival`].
#' @param ... Further arguments passed to or from other methods.
#' @return 
#' A `stateprobs` object.
#' @examples 
#' library("data.table")
#' library("survival")
#' 
#' # This example shows how to simulate a partitioned survival model by
#' # manually constructing a "survival" object. We will consider a case in which
#' # Cox proportional hazards models from the survival package---which are not
#' # integrated with hesim---are used for parameter estimation. We will use 
#' # point estimates in the example, but bootstrapping, Bayesian modeling,
#' # or other techniques could be used to draw samples for a probabilistic 
#' # sensitivity analysis. 
#' 
#' # (0) We first setup our model per usual by defining the treatment strategies,
#' # target population, and health states
#' hesim_dat <- hesim_data(
#'   strategies = data.table(strategy_id = 1:3,
#'                           strategy_name = c("SOC", "New 1", "New 2")),
#'   patients = data.table(patient_id = 1:2,
#'                         female = c(0, 1),
#'                         grp_id = 1),
#'   states = data.table(state_id = 1:2,
#'                       state_name = c("Stable", "Progression"))
#' )
#'
#' # (1) Next we will estimate Cox models with survival::coxph(). We illustrate 
#' # by predicting progression free survival (PFS) and overall survival (OS)
#' ## Fit models
#' onc3_pfs_os <- as_pfs_os(onc3, patient_vars = c("patient_id", "female",
#'                                                 "strategy_name"))
#' fit_pfs <- coxph(Surv(pfs_time, pfs_status) ~ strategy_name + female,
#'                  data = onc3_pfs_os)
#' fit_os <- coxph(Surv(os_time, pfs_status) ~ strategy_name + female,
#'                 data = onc3_pfs_os)
#'
#' ## Predict survival on input data
#' surv_input_data <- expand(hesim_dat)
#' times <- seq(0, 14, 1/12)
#' predict_survival <- function(object, newdata, times) {
#'   surv <- summary(survfit(object, newdata = newdata, se.fit = FALSE),
#'                   t = times)
#'   pred <- newdata[rep(seq_len(nrow(newdata)), each = length(times)), ]
#'   pred[, sample := 1] # Point estimates only in this example
#'   pred[, time := rep(surv$time, times = nrow(newdata))]
#'   pred[, survival := c(surv$surv)]
#'   return(pred[, ])
#' }
#' pfs <- predict_survival(fit_pfs, newdata = surv_input_data, times = times)
#' os <- predict_survival(fit_os, newdata = surv_input_data, times = times)
#' surv <- rbind(
#'   as.data.table(pfs)[, curve := 1L],
#'   as.data.table(os)[, curve := 2L]
#' )
#'
#' ## Convert predictions to a survival object
#' surv <- survival(surv, t = "time")
#' \dontrun{autoplot(surv)}
#' 
#' # (2) We can then compute state probabilities from the survival object
#' stprobs <- sim_stateprobs(surv)
#' 
#' # (3) Finally, we can use the state probabilities to compute QALYs and costs
#' ## A dummy utility model to illustrate
#' utility_tbl <- stateval_tbl(
#'   data.table(state_id = 1:2,
#'              est = c(1, 1)
#'   ),
#'   dist = "fixed"
#' )
#' utilitymod <- create_StateVals(utility_tbl, 
#'                                hesim_data = hesim_dat,
#'                                n = 1)
#'
#' ## Instantiate Psm class and compute QALYs
#' psm <- Psm$new(utility_model = utilitymod)
#' psm$stateprobs_ <- stprobs
#' psm$sim_qalys()
#' psm$qalys_
#' 
#' @seealso [`survival`]
#' @export
sim_stateprobs.survival <- function(x, ...) {
  state_id <- NULL
  
  # Size and attributes
  if (is.null(attr(x, "size"))) stop("'size' attribute missing from 'x'.")
  n_samples <- attr(x, "size")[["n_samples"]]
  n_strategies <- attr(x, "size")[["n_strategies"]]
  n_patients <- attr(x, "size")[["n_patients"]]
  n_states <- attr(x, "size")[["n_states"]]
  n_curves <- n_states - 1
  n_times <- attr(x, "size")[["n_times"]]
  N <- n_samples * n_strategies * n_patients
  unique_times <- x$t[1:n_times]
  
  
  # Simulate state probabilities
  ## Reshape into array for easier handling with C++
  x2 <- array(x$survival, dim = c(n_times,  n_curves, N))
  
  ## Compute state probabilities via partitioned survival analysis
  stprobs <- C_psm_sim_stateprobs(x2)
  
  # Link simulated state probabilities to ID variables
  ## First expand 'x' to include rows for each health state 
  ## (n_states = n_curves + 1)
  rep_every_n <- function(i, n, times = 2) {
    v <- rep(1, length(i))
    v[seq_along(v) %% n == 0L] <- times
    rep(i, v)
  }
  n_curves * n_times
  new_rows <- rep_every_n(1:nrow(x), n = n_curves * n_times, times = n_times + 1)
  cols <- c("sample", "strategy_id", "patient_id", "grp_id", "patient_wt")
  cols <- cols[cols %in% colnames(x)]
  stprobs_df <- x[new_rows, cols, with = FALSE]
  stprobs_df[, state_id := rep(rep(1:n_states, each = n_times),
                               times = N)]
  stprobs_df[, t := rep(unique_times, times = N * n_states)]
  
  ## Then add the probabilities as a column
  stprobs_df[, ("prob") := c(stprobs$prob)]
  
  # Summarize curve crossing
  nr <- nrow(stprobs_df)
  n_cross <- sum(stprobs$cross)
  if (n_cross > 0) {
    percent_cross <- paste0(formatC(n_cross/nr * 100, format = "f", digits = 1), 
                            "%")
    warning(paste0("The survival curves were crossed ",
                   n_cross, "/", nr, " (", percent_cross, ") ",
                   "of the time."),
            call. = FALSE)
  }
  
  # Create object and set attributes
  setattr(stprobs_df, "class", 
          c("stateprobs", "data.table", "data.frame"))
  setattr(stprobs_df, "size", 
          c(n_samples = n_samples, n_strategies = n_strategies, 
            n_patients = n_patients, n_states = n_states, n_times = n_times))
  return(stprobs_df[, ])
}

# Cost and QALYs ---------------------------------------------------------------
#' Costs object
#'
#' An object of class `costs` returned from methods 
#' `$sim_costs()` in model classes that store simulated costs. 
#' 
#' @section Components:
#' A `costs` object inherits from `data.table` and contains
#' the following columns:
#' 
#' \describe{
#'   \item{sample}{A random sample from the PSA.}
#'   \item{strategy_id}{The treatment strategy ID.}
#'   \item{patient_id}{The patient ID.}
#'   \item{state_id}{The health state ID.}
#'   \item{dr}{The rate used to discount costs.}
#'   \item{category}{The cost category (e.g., drug costs, medical costs, etc).}
#'   \item{costs}{The simulated cost values.}
#' }
#'
#' @name costs
NULL

#' Quality-adjusted life-years object
#'
#' An object of class `qalys` returned from methods 
#' `$sim_qalys()` in model classes that store simulated 
#' quality-adjusted life-years (QALYs).
#' 
#' @section Components:
#' A `qalys` object inherits from `data.table` and contains
#' the following columns:
#' 
#' \describe{
#'   \item{sample}{A random sample from the PSA.}
#'   \item{strategy_id}{The treatment strategy ID.}
#'   \item{patient_id}{The patient ID.}
#'   \item{state_id}{The health state ID.}
#'   \item{dr}{The rate used to discount QALYs.}
#'   \item{category}{A single category always equal to "qalys".}
#'   \item{qalys}{The simulated values of QALYs.}
#' }
#' If the argument `lys = TRUE`, then the `data.table` also contains a column
#' `lys` containing simulated life-years.
#' @name qalys
NULL

sim_ev <- function (object, ...) {
  UseMethod("sim_ev", object)
}

sim_ev.NULL <- function(object, ...) {
  if (is.null(object)) {
    stop("You must first simulate state probabilities using '$sim_stateprobs'.",
         call. = FALSE)
  }
}


sim_ev.stateprobs <- function(object, statevalmods, categories, dr = .03,
                              integrate_method = c("trapz", "riemann_left", "riemann_right")){
  integrate_method <- match.arg(integrate_method)
  state_id <- NULL
  
  # Checks
  ## State probabilities
  if(is.null(object)){
    stop("You must first simulate health state probabilities.",
         call. = FALSE)
  }
  
  ## Discount rate
  check_dr(dr)
  
  ## The size of ID variables is correct
  check_StateVals(statevalmods, object)
  
  # Simulate
  res <- data.table(C_sim_ev(object[state_id != max(state_id)],
                             statevalmods,
                             dr, categories,
                             unique(object$t),
                             integrate_method))
  res[, sample := sample + 1]
  if (!"patient_wt" %in% colnames(object)) res[, ("patient_wt") := NULL]
  return(res[])
} 

sim_los <- function (object, utility_model, dr, 
                     integrate_method = c("trapz", "riemann_left", "riemann_right")) {
  state_id <- NULL # To avoid no visible binding note
  n_samples <- utility_model$params$n_samples
  id <- get_id_object(utility_model)
  n_obs <- n_samples * id$n_strategies * id$n_patients * id$n_states
  los <- C_sim_los(stateprobs = object[state_id != max(state_id)]$prob,
                   n_obs = n_obs,
                   dr = dr,
                   times = unique(object$t),
                   integrate_method = match.arg(integrate_method))
  return(los)
}

#' Expected values
#' 
#' Simulate costs and quality-adjusted life-years (QALYs) as a function of
#' simulated state occupancy probabilities. 
#' 
#' @param object A [`stateprobs`] object.
#' @param utility_model A single object of class [`StateVals`] used
#' to simulate utility.
#' @param cost_models A list of objects of class [`StateVals`] used
#' to simulate costs.
#' @param dr Discount rate. 
#' @param integrate_method Method used to integrate state values when computing 
#' weighted length of stay. Options are `trapz` for the trapezoid rule,
#' `riemann_left` left for a left Riemann sum, and  
#' `riemann_right` right for a right Riemann sum.
#' @param lys If `TRUE`, then life-years are simulated in addition to 
#' QALYs. 
#' @keywords internal
#' @return [`sim_costs()`] and [`sim_qalys()`] return objects of class
#' [`costs`] and [`qalys`], respectively. 
#' @details 
#' See `vignette("expected-values")` for details.
#'
#' @name sim_ev
sim_qalys <- function(object, utility_model, dr, integrate_method, lys){
  utility_model$check()
  qalys <- sim_ev(object,
                  list(utility_model),
                  "qalys",
                  dr,
                  integrate_method)
  if (lys){
    los <- sim_los(object,
                   utility_model,
                   dr,
                   integrate_method)
    qalys[, lys := los]
  }
  qalys[, ("category") := NULL]
  setnames(qalys, "value", "qalys")
  setattr(qalys, "class", 
          c("qalys", "data.table", "data.frame"))
  return(qalys[, ])
}

#' @rdname sim_ev
sim_costs <- function(object, cost_models, dr, integrate_method){
  if(!is.list(cost_models)){
    stop("'cost_models' must be a list", call. = FALSE)
  }
  for (i in 1:length(cost_models)){
    cost_models[[i]]$check()
  }
  if (is.null(names(cost_models))){
    categories <- paste0("Category ", seq(1, length(cost_models)))
  } else{
    categories <- names(cost_models)
  }   
  costs <- sim_ev(object,
                  cost_models,
                  categories,
                  dr,
                  integrate_method)
  setnames(costs, "value", "costs")
  setattr(costs, "class", 
          c("costs", "data.table", "data.frame"))
  return(costs[, ])
}

# Cost-effectiveness object ----------------------------------------------------
#'  A cost-effectiveness object
#'
#' An object that summarizes simulated measures of clinical effectiveness and 
#' costs from a simulation model for use in a cost-effectiveness analysis.
#'
#' 
#' @format 
#' A list containing two elements:
#' \itemize{
#' \item{`costs`}{ Total (discounted) costs by category.}
#' \item{`qalys`}{ (Discounted) quality-adjusted life-years.}
#' }
#' 
#' @section Costs:
#' The `costs` `data.table` contains the following columns:
#' \describe{
#' \item{category}{The cost category.}
#' \item{dr}{The discount rate.}
#' \item{sample}{A randomly sampled parameter set from the probabilistic sensitivity analysis (PSA)}
#' \item{strategy_id}{The treatment strategy ID.}
#' \item{grp_id}{An optional column denoting a subgroup. If not included, it is
#'  assumed that a single subgroup is being analyzed.}
#' \item{costs}{Costs.}
#' }
#' 
#' @section Quality-adjusted life-years:
#' The `qalys` `data.table` contains the following columns:
#' \describe{
#' \item{dr}{The discount rate.}
#' \item{sample}{A randomly sampled parameter set from the probabilistic sensitivity analysis (PSA)}
#' \item{strategy_id}{The treatment strategy ID.}
#' \item{grp_id}{An optional column denoting a subgroup. If not included, it is 
#' assumed that a single subgroup is being analyzed.}
#' \item{qalys}{Quality-adjusted life-years}
#' }
#' 
#' @name ce
NULL