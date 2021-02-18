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
#' An object of class `survival` returned from `Psm$sim_survival()` or 
#' `PsmCurves$survival()`.
#' 
#' @section Components:
#' A `survival` object inherits from `data.table` and contains
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
#'
#' @seealso [`Psm`], [`PsmCurves`]    
#' @name survival
NULL

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