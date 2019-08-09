# transprob_tbl ----------------------------------------------------------------
#' Table to store transition probabilities 
#' 
#' Create a table for storing transition probabilities for a discrete time state
#' transition model stratified by treatment strategy, patient, transition, and
#' (optionally) time interval.
#' 
#' @param tbl A \code{data.frame} or \code{data.table} for storing parameter 
#' values. See "Details" for specifics. 
#' @param dist Probability distribution used to sample parameters for a 
#' probabilistic sensitivity analysis (PSA). 
#' @param hesim_data A \code{\link{hesim_data}} object. Required to specify 
#' patients if not included as columns in \code{tbl}. Not required if \code{tbl} 
#' includes one row for each treatment strategy, patient, and transition combination.
#' @details 
#' \code{tbl} is a \code{data.table} that must contain columns for
#' treatment strategies (\code{strategy_id}) and health state 
#' transitions (\code{transition_id}). The number of total transitions
#' must be equal to the square of number of health states. The transition
#' IDs are, by default, integers starting at 1 and ordered rowwise in 
#' transition matrices where each (i,j) pair denotes a transition from
#' state i to state j. For example, in a model with two health states 
#' transition 1 would be the transition from state 1 to itself, transition 2 would
#' be the transition from state 1 to state 2, transition 3 would be the
#' transition from state 2 to state 1, and transition 4 would be the transition 
#' from state 2 to itself. Alternatively, transition IDs can be specified in
#' \code{hesim_data} by including a \code{transitions} table. Each row is a 
#' treatment strategy and transition pair. 
#' 
#' Transition probabilities may also vary by patient (\code{patient_id}) 
#' and time intervals (\code{time_start}). If both of these variables are included, 
#' then each row is a treatment strategy, transition, patient, and time interval pair.
#' If no patient ID is specified in the table, then the patients in the simulation
#' can be specified using the \code{hesim_data} argument and it will be assumed
#' that transition probabilities are the same across patients. 
#' 
#' \code{tbl} must also contain columns summarizing the state values for each
#' row. Parameter values for a PSA can be drawn from Dirichlet distributions 
#' (separate for each treatment strategy and optionally patient subgroup 
#' and time interval) with \code{dist = 'dirichlet'}. Alternatively, 
#' \code{fixed} can be used if estimates are known with certainty
#' and \code{custom} can be used if transition probabilities for a PSA
#'  have been previously sampled from an arbitrary probability distribution.
#'  If \code{dist = "custom"}, then \code{tbl} must include a column 
#'  named \code{sample} (an integer vector denoting a unique random draw).
#'  The columns in \code{tbl} related to the values of the 
#'  transition probabilities depend on the value of \code{dist}. They are
#'  as follows:
#' 
#' \describe{
#' \item{fixed}{\code{est}}
#' \item{dirichlet}{\code{alpha}}
#' \item{custom}{\code{value}}
#' }
#'  
#' @return An object of class "transprob_tbl", which is a \code{data.table} of
#' parameter values with attributes for \code{dist} and optionally 
#'  \code{patients}. \code{tbl} is in the same format as described in "Details". 
#'  \code{patients} is a 
#' \code{data.table} with at leat one column containing \code{patient_id}.
#' @seealso \code{\link{create_CohortDtstmTrans}},  \code{\link{CohortDtstmTrans}}
#' @export
transprob_tbl <- function(tbl, 
                          dist = c("custom", "dirichlet", "fixed"),
                          hesim_data = NULL){
  x <- ParamsTbl$new(data.table(tbl), match.arg(dist), hesim_data)
  x$create_tbl("transprob_tbl")
  return(x$tbl[, ])  
}

# CohortDtstmTrans -------------------------------------------------------------
#' @export
CohortDtstmTrans <- R6::R6Class("CohortCtstmTrans",
  public = list(
    input_mats = NULL,
    params = NULL,
    start_stateprobs = NULL,
    
    initialize = function(input_mats, params, start_stateprobs){
      self$input_mats <- input_mats
      self$params <- params
      self$start_stateprobs <- start_stateprobs
    },
    
    sim_stateprobs = function(n_cycles){
      return(2)
    }
      
  )
)

# create_CohortDtstmTrans ------------------------------------------------------
#' Create \code{CohortDtstmTrans} object
#' 
#' \code{create_CohortDtstmTrans} creates an object of class
#'  \code{\link{CohortDtstmTrans}} from a \code{\link{transprob_tbl}}
#'  object. 
#' @param object An object of class "transprob_tbl".
#' @param n Number of random observations of the parameters to draw for the PSA 
#' if the parameters have not previously been sampled. See "details". 
#' @param start_stateprobs Initial state probabilities. By default all patients
#' are assumed to be in the first health state.
#' @param ... Further arguments passed to or from other methods. Currently unused. 
#' @details The argument \code{n} is not used if the option \code{dist = "custom"} 
#'  is used in \code{\link{transprob_tbl}} since the parameters have already be 
#'  drawn for the PSA.
#' @return Returns an \code{\link{R6Class}} object of class \code{\link{CohortDtstmTrans}}.
#' @seealso \code{\link{CohortDtstmTrans}}
#' @export
create_CohortDtstmTrans <- function(object, ...){
  UseMethod("create_CohortDtstmTrans", object)
} 

#' @rdname create_CohortDtstmTrans
#' @export  
create_CohortDtstmTrans.transprob_tbl <- function(object, n = 1000, 
                                                  start_stateprobs = NULL, ...){
  x <- CreateFromParamsTbl$new(object, n)
  x$prep()
  if (is.null(start_stateprobs)){
    start_stateprobs <- c(1, rep(0, x$n_states - 1))
  } 
  transmod <- CohortDtstmTrans$new(input_mats = x$input_mats, params = x$params,
                                   start_stateprobs = start_stateprobs)
  return(transmod)
}