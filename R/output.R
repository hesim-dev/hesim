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

# Survival object --------------------------------------------------------------
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

# State probability object -----------------------------------------------------
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

# Costs object -----------------------------------------------------------------
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

# QALYs object -----------------------------------------------------------------
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