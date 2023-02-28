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
#' The object also contains `size` and `absorbing` attributes.
#' The `size` attribute is a numeric vector with the elements `n_samples`,
#' `n_strategies`, `n_patients`, and `n_states` denoting the number of samples,
#' treatment strategies, patients, and health states The `absorbing` attribute
#' is a numeric vector containing the absorbing health states; i.e., the
#' health states that cannot be transitioned from. Operationally, an
#' absorbing state is a row in a transition matrix (as in the `trans_mat` field
#' of the `IndivCtstmTrans` class) with all `NA`s.
#'
#' @seealso A disease progression object can be simulated with either the
#' [`IndivCtstm`] or [`IndivCtstmTrans`] classes.
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
  if (prod(id_n) != nrow(x)) {
    stop(paste0(
      "The number of rows in 'data' must be equal to the product of the ",
      "number of unique values of the 'sample', 'strategy_id', 'patient_id' ",
      "'grp_id', 'curve', and 't' columns."
    ))
  }

  # Return
  setattr(x, "class", c("survival", "data.table", "data.frame"))
  setattr(
    x, "size",
    c(
      n_samples = id_n[["sample"]],
      n_strategies = id_n[["strategy_id"]],
      n_patients = id_n[["patient_id"]],
      n_states = id_n[["curve"]] + 1,
      n_times = id_n[["t"]]
    )
  )
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
#' In cohort models, the object also contains `size` and `absorbing` attributes.
#' The `size` attribute is a numeric vector with the elements `n_samples`,
#' `n_strategies`, `n_patients`, `n_states`, and
#' `n_times` denoting the number of samples, treatment strategies, patients,
#'  health states, and times. The `absorbing` attribute is a numeric vector
#'  containing the absorbing health states (see the `absorbing` field of the
#'  [`CohortDtstmTrans`] class for more details).
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
#' @keywords internal
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
#' @details In an \eqn{N}-state partitioned survival model there are \eqn{N-1} survival curves
#' and \eqn{S_n(t)} is the cumulative survival function denoting the probability of
#' survival to health state \eqn{n} or a lower indexed state beyond time \eqn{t}.
#' The probability that a patient is in health state 1 at time \eqn{t} is simply
#' \eqn{S_1(t)}. State membership in health states \eqn{2,\ldots, N -1} is calculated
#' as \eqn{S_n(t) - S_{n-1}(t)}. Finally, the probability of being in the final
#' health state \eqn{N} (i.e., the death state) is \eqn{1-S_{N-1}(t)}, or
#' one minus the overall survival curve.
#'
#' In some cases, the survival curves may cross. `hesim` will issue a warning
#' but the function will still run. Probabilities will be set to 0 in a health state
#' if the prior survival curve lies above the curve for state \eqn{n};
#' that is, if \eqn{S_n(t) < S_{n-1}(t)}, then the probability of being in state \eqn{n}
#' is set to 0 and \eqn{S_n(t)} is adjusted to equal \eqn{S_{n-1}(t)}. The
#' probability of being in the final health state is also adjusted if necessary to
#' ensure that probabilities sum to 1.
#'
#' @return
#' A [`stateprobs`] object.
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
#'   strategies = data.table(
#'     strategy_id = 1:3,
#'     strategy_name = c("SOC", "New 1", "New 2")
#'   ),
#'   patients = data.table(
#'     patient_id = 1:2,
#'     female = c(0, 1),
#'     grp_id = 1
#'   ),
#'   states = data.table(
#'     state_id = 1:2,
#'     state_name = c("Stable", "Progression")
#'   )
#' )
#'
#' # (1) Next we will estimate Cox models with survival::coxph(). We illustrate
#' # by predicting progression free survival (PFS) and overall survival (OS)
#' ## Fit models
#' onc3_pfs_os <- as_pfs_os(onc3, patient_vars = c(
#'   "patient_id", "female",
#'   "strategy_name"
#' ))
#' fit_pfs <- coxph(Surv(pfs_time, pfs_status) ~ strategy_name + female,
#'   data = onc3_pfs_os
#' )
#' fit_os <- coxph(Surv(os_time, pfs_status) ~ strategy_name + female,
#'   data = onc3_pfs_os
#' )
#'
#' ## Predict survival on input data
#' surv_input_data <- expand(hesim_dat)
#' times <- seq(0, 14, 1 / 12)
#' predict_survival <- function(object, newdata, times) {
#'   surv <- summary(survfit(object, newdata = newdata, se.fit = FALSE),
#'     t = times
#'   )
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
#' \dontrun{
#' autoplot(surv)
#' }
#'
#' # (2) We can then compute state probabilities from the survival object
#' stprobs <- sim_stateprobs(surv)
#'
#' # (3) Finally, we can use the state probabilities to compute QALYs and costs
#' ## A dummy utility model to illustrate
#' utility_tbl <- stateval_tbl(
#'   data.table(
#'     state_id = 1:2,
#'     est = c(1, 1)
#'   ),
#'   dist = "fixed"
#' )
#' utilitymod <- create_StateVals(utility_tbl,
#'   hesim_data = hesim_dat,
#'   n = 1
#' )
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
  x2 <- array(x$survival, dim = c(n_times, n_curves, N))

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
    times = N
  )]
  stprobs_df[, t := rep(unique_times, times = N * n_states)]

  ## Then add the probabilities as a column
  stprobs_df[, ("prob") := c(stprobs$prob)]

  # Summarize curve crossing
  nr <- nrow(stprobs_df)
  n_cross <- sum(stprobs$cross)
  if (n_cross > 0) {
    percent_cross <- paste0(
      formatC(n_cross / nr * 100, format = "f", digits = 1),
      "%"
    )
    warning(paste0(
      "The survival curves were crossed ",
      n_cross, "/", nr, " (", percent_cross, ") ",
      "of the time."
    ),
    call. = FALSE
    )
  }

  # Create object and set attributes
  setattr(
    stprobs_df, "class",
    c("stateprobs", "data.table", "data.frame")
  )
  setattr(
    stprobs_df, "size",
    c(
      n_samples = n_samples, n_strategies = n_strategies,
      n_patients = n_patients, n_states = n_states, n_times = n_times
    )
  )
  setattr(stprobs_df, "absorbing", n_states)
  return(stprobs_df[, ])
}

# Cost and QALY objects --------------------------------------------------------
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
#'   \item{grp_id}{The subgroup ID.}
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
#'   \item{grp_id}{The subgroup ID.}
#'   \item{state_id}{The health state ID.}
#'   \item{dr}{The rate used to discount QALYs.}
#'   \item{category}{A single category always equal to "qalys".}
#'   \item{qalys}{The simulated values of QALYs.}
#' }
#' If the argument `lys = TRUE`, then the `data.table` also contains a column
#' `lys` containing simulated life-years.
#' @name qalys
NULL

# Simulate expected values (inclusive of costs and QALYs) ----------------------
#' @export
sim_ev <- function(object, ...) {
  UseMethod("sim_ev", object)
}

sim_ev.NULL <- function(object, ...) {
  if (is.null(object)) {
    stop("You must first simulate state probabilities using '$sim_stateprobs'.",
      call. = FALSE
    )
  }
}

#' Expected values from state probabilities
#'
#' Simulate expected values as a function of simulated state occupancy
#' probabilities, with simulation of costs and quality-adjusted life-years
#' (QALYs) as particular use cases.
#'
#' @param object A [`stateprobs`] object.
#' @param model,models An object or list of objects of class [`StateVals`] used to
#' model state values. When using `sim_qalys()`, this should be
#' a single model for utility. With `sim_costs()`, a list of models should be
#' used with one model for each cost category. Finally, with `sim_ev()`,
#' this may either be a single model or a list of models. May also be `NULL`,
#' in which case length of stay is computed based on the state probabilities
#' contained in `object`.
#' @param dr Discount rate.
#' @param integrate_method Method used to integrate state values when computing
#' costs or QALYs. Options are `trapz` (the default) for the trapezoid rule,
#' `riemann_left` left for a left Riemann sum, and
#' `riemann_right` right for a right Riemann sum.
#' @param outcome_name Name of the column indicating the outcome corresponding
#' to each model. Only used if `models` is a list. Default is `"outcome"`.
#' @param value_name Name of the column containing values of the outcome. Default
#' is `"value"`.
#' @param lys If `TRUE`, then life-years are simulated in addition to
#' QALYs.
#' @param ... Currently unused.
#'
#' @details
#' Expected values in cohort models (i.e.,  those implemented with
#' the [`CohortDtstm`] and [`Psm`] classes) are mean outcomes for patients comprising
#' the cohort. The method used to simulate expected values depends on the
#' `$method` field in the [`StateVals`] object(s) stored in `model(s)`. If
#' `$method = "starting"`, then state values represent a one-time value that
#' occurs at time 0.
#'
#' The more common use case is `$method = "wlos"`, or a "weighted length of stay".
#' That is, expected values for each health state can be thought of as state values
#' weighted by the time a patient spends in each state (and discounted by a
#' discount factor that depends on the discount rate `dr`). The
#' precise computation proceeds in four steps. In the first step, the probability
#' of being in each health state at each discrete time point is simulated
#' (this is the output contained in the [`stateprobs`] object). Second, a
#' [`StateVals`] model is used to predict state values at each time point.
#' Third an expected value at each time point is computed by multiplying the
#' state probability, the state value, and the discount factor. Fourth, the
#' expected values at each time point are summed across all time points.
#'
#' The summation in the fourth step can be thought of as a discrete approximation
#' of an integral. In particular, the limits of integration can be partitioned
#' into time intervals, with each interval containing a start and an end.
#' The `integrate_method` argument determines the approach used
#' for this approximation:
#'
#' 1. A left Riemann sum (`integrate_method = "riemann_left"`) uses expected values
#' at the start of each time interval.
#' 2. A right Riemann sum (`integrate_method = "riemann_right"`) uses expected values
#' at the end of each time interval.
#' 3. The trapezoid rule (`integrate_method = "trapz"`) averages expected values
#' at the start and end of each time interval. (This will generally be the
#'  most accurate and is recommended.)
#'
#' Mathematical details are provided in the reference within the "References"
#' section below.
#'
#' @note The ID variables in the state value models in `models` must be
#' consistent with the ID variables contained in `object`. In particular,
#' the `models` should predict state values for each non-absorbing health state
#' in `object`; that is, the number of health states modeled with the
#'  `models` should equal the number of health states in `object` less the number
#'  of absorbing states.
#'
#' The absorbing states are saved as an attribute named `absorbing` to
#' [`stateprobs`] objects. When simulating state probabilities with a
#' [`CohortDtstmTrans`] object, the absorbing state is determined by the
#' `absorbing` field in the class; in a `Psm` (or with
#' [sim_stateprobs.survival()]), the absorbing state is always equal to the
#' final health state.
#'
#' @return `sim_ev()` returns a `data.table` with the following columns:
#' \describe{
#'   \item{sample}{A random sample from the PSA.}
#'   \item{strategy_id}{The treatment strategy ID.}
#'   \item{patient_id}{The patient ID.}
#'   \item{grp_id}{The subgroup ID.}
#'   \item{state_id}{The health state ID.}
#'   \item{dr}{The rate used to discount costs.}
#'   \item{outcome}{The outcome corresponding to each model in `models`.
#'   Only included if `models` is a list.}
#'   \item{value}{The expected value.}
#' }
#' The names of the `outcome` and `value` columns may be changed with the
#' `value_name` and `outcome_name` arguments. `sim_costs()` and `sim_qalys()`
#' return similar objects, that are of class [`costs`] and [`qalys`], respectively.
#'
#' @seealso State probabilities can be simulated using the
#' `$sim_stateprobs()` methods from either the [`CohortDtstmTrans`]
#' (or [`CohortDtstm`]) or [`Psm`] classes. State probabilities can also be
#' computed directly from survival curves with the generic method
#' [sim_stateprobs.survival()].
#'
#' Costs and QALYs are typically computed within the `R6` model classes
#' using the `$sim_costs()` and `$sim_qalys()` methods. For instance, see the
#' documentation and examples for the [`CohortDtstm`] and [`Psm`] classes.
#' The `sim_qalys()` and `sim_costs()` functions are exported to give users
#' additional flexibility when creating their own modeling pipelines.
#' `sim_ev()` may be useful for computing outcomes other than costs or QALYs.
#'
#' [`costs`] and [`qalys`] objects can be passed to [summarize_ce()] to
#' create a cost-effectiveness object for performing a cost-effectiveness analysis
#' with [cea()]. Although note that typically the `$summarize()` method
#' belonging to the [`CohortDtstm`] or [`Psm`] classes would be used instead.
#'
#' Use the [`IndivCtstm`] class to simulate costs and QALYs with an individual
#' continuous-time state transition model.
#'
#' @example man-roxygen/example-sim_ev.stateprobs.R
#'
#' @references [Incerti and Jansen (2021)](https://arxiv.org/abs/2102.09437).
#' See Section 2.1 for mathematical details.
#'
#' @export
#' @name sim_ev
#' @aliases sim_costs sim_qalys
sim_ev.stateprobs <- function(object, models = NULL, dr = .03,
                              integrate_method = c("trapz", "riemann_left", "riemann_right"),
                              value_name = "value", outcome_name = "outcome",
                              ...) {
  integrate_method <- match.arg(integrate_method)
  state_id <- NULL

  # Some standard checks
  if (is.null(object)) {
    stop("You must first simulate health state probabilities.",
      call. = FALSE
    )
  }
  check_dr(dr)

  # Resize the state probabilities based on the absorbing states
  absorbing <- attr(object, "absorbing")
  if (is.null(absorbing)) {
    stprobs <- object
  } else {
    stprobs <- object[!state_id %in% absorbing]
  }

  # Case where models is NULL
  if (is.null(models)) {
    out <- data.table(
      stprobs[t == 0][, c("prob", "t") := NULL],
      value = sim_los(
        object,
        dr,
        integrate_method
      )
    )
    setnames(out, "value", value_name)
    return(out)
  }

  # Case where models is not NULL
  ## Differentiate between a single model or a list of models being provided
  if (inherits(models, "StateVals")) {
    models <- list(models)
    model_is_list <- FALSE
  } else {
    model_is_list <- TRUE
  }

  ## Add default names if required
  if (is.null(names(models))) {
    names(models) <- paste0("Outcome ", seq(1, length(models)))
  }

  ## Check the state value models, particularly that the size of the
  # ID variables is correct
  check_state_vals(models, object)


  ## Simulate
  out <- data.table(C_sim_ev(
    stprobs,
    models,
    dr,
    names(models),
    unique(object$t),
    integrate_method
  ))
  if (!model_is_list) {
    out[, ("outcome") := NULL]
  } else {
    setnames(out, "outcome", outcome_name)
  }
  out[, sample := sample + 1]
  if (!"patient_wt" %in% colnames(object)) out[, ("patient_wt") := NULL]
  setnames(out, "value", value_name)
  return(out[])
}

sim_los <- function(object, dr,
                    integrate_method = c("trapz", "riemann_left", "riemann_right")) {
  state_id <- NULL # To avoid no visible binding note
  absorbing <- attr(object, "absorbing")
  if (is.null(absorbing)) {
    stprobs <- object
  } else {
    stprobs <- object[!state_id %in% absorbing]
  }
  unique_times <- unique(object$t)
  los <- C_sim_los(
    stateprobs = stprobs$prob,
    n_obs = nrow(stprobs) / length(unique_times),
    dr = dr,
    times = unique_times,
    integrate_method = match.arg(integrate_method)
  )
  return(los)
}

#' @rdname sim_ev
#' @export
sim_qalys <- function(object, model, dr = .03,
                      integrate_method = c("trapz", "riemann_left", "riemann_right"),
                      lys = TRUE) {
  model$check()
  qalys <- sim_ev(object,
    models = model,
    dr = dr,
    integrate_method = integrate_method,
    value_name = "qalys"
  )
  if (lys) {
    los <- sim_los(
      object,
      dr,
      integrate_method
    )
    qalys[, lys := los]
  }
  setattr(
    qalys, "class",
    c("qalys", "data.table", "data.frame")
  )
  return(qalys[, ])
}

#' @rdname sim_ev
#' @export
sim_costs <- function(object, models, dr = .03,
                      integrate_method = c("trapz", "riemann_left", "riemann_right")) {
  if (!is.list(models)) {
    stop("'models' must be a list", call. = FALSE)
  }
  if (is.null(names(models))) {
    names(models) <- paste0("Category ", seq(1, length(models)))
  }
  costs <- sim_ev(object,
    models = models,
    dr = dr,
    integrate_method = integrate_method,
    value_name = "costs",
    outcome_name = "category"
  )
  setattr(
    costs, "class",
    c("costs", "data.table", "data.frame")
  )
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
