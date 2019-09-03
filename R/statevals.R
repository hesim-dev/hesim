# stateval_tbl -----------------------------------------------------------------
#' Table to store state value parameters
#' 
#' Create a table for storing parameter estimates used to simulate costs or 
#' utility in an economic model by treatment strategy, patient, health state, and
#' (optionally) time interval. 
#' 
#' @param tbl A \code{data.frame} or \code{data.table} for storing parameter 
#' values. See "Details" for specifics. 
#' @param dist Probability distribution used to sample parameters for a 
#' probabilistic sensitivity analysis (PSA). 
#' @param hesim_data A \code{\link{hesim_data}} object. Required to specify 
#' treatment strategies , patients,
#'  and/or health states not included as columns
#' in \code{tbl}, or, to match patients in \code{tbl} to groups.
#'  Not required if \code{tbl} includes one row for each treatment strategy, patient, and
#' health state combination. Patients are matched to groups by specifying both a 
#' \code{patient_id} and a \code{grp_id} column in the \code{patients} table;
#' if \code{hesim_data = NULL} but \code{grp_id} is included as a column
#' in \code{tbl}, then each group is assumed to be a unique patient. 
#' 
#' @details 
#' \code{tbl} is a \code{data.table} containing columns for treatment 
#' strategies (\code{strategy_id}), patient subgroups (\code{grp_id}),
#' health states (\code{state_id}), and/or the start of time intervals 
#' (\code{time_start}). The table must contain at least one column
#' named \code{strategy_id}, \code{grp_id}, or \code{state_id}, 
#' but does not need to contain all of them. Each row denotes a unique 
#' treatment strategy, patient subgroup, health state, and/or time interval pair.
#' 
#' 
#' \code{tbl} must also contain columns summarizing the state values for each
#' row, which depend on the probability distribution selected with \code{dist}. 
#' Available distributions include the normal (\code{norm}), beta (\code{beta}),
#' gamma (\code{gamma}), lognormal (\code{lnorm}), and uniform (\code{unif})
#'  distributions. In addition, the option \code{fixed} can be used if estimates
#'  are known with certainty and \code{custom} can be used if 
#'  parameter values for a PSA  have been previously
#' sampled from an arbitrary probability distribution.
#'  The columns in \code{tbl} that must be included,
#'  by distribution, are:
#' 
#' \describe{
#' \item{norm}{\code{mean} and \code{sd}}
#' \item{beta}{\code{mean} and \code{se} or \code{shape1} and \code{shape2}}
#' \item{gamma}{\code{mean} and \code{se}, \code{shape} and \code{rate}, 
#' or \code{shape} and {scale}}
#' \item{lnorm}{\code{meanlog} or \code{sdlog}}
#' \item{unif}{\code{min} and \code{max}}
#' \item{fixed}{\code{est}}
#' \item{custom}{\code{sample} and \code{value}}
#' }
#' 
#' Note that if \code{dist = "custom"}, then \code{tbl} must include a column 
#' named \code{sample} (an integer vector denoting a unique random draw) and
#'  \code{value} (denoting the value of the randomly sampled parameter). In this case, there is a unique
#' row in \code{tbl} for each random draw (\code{sample}) and
#'  each combination of strategies, patients, health states, and/or time intervals.
#' Again, \code{tbl} must contain at least one column
#' named \code{strategy_id}, \code{grp_id}, or \code{state_id},
#'  but does not need to contain them all.
#'  
#'  
#' 
#' @return An object of class "stateval_tbl", which is a \code{data.table} of
#' parameter values with attributes for \code{dist} and optionally 
#' \code{strategy_id}, \code{patients}, and \code{state_id}. \code{tbl} 
#' is in the same format as described in "Details". \code{patients} is a 
#' \code{data.table} with one column containing \code{patient_id} and 
#' optionally a second column containing \code{grp_id}.
#' @seealso \code{\link{create_StateVals}}, \code{\link{StateVals}}
#' @examples 
#' strategies <- data.frame(strategy_id = c(1, 2))
#' patients <- data.frame(patient_id = seq(1, 3),
#'                        grp_id = c(1, 1, 2),
#'                        age = c(45, 50, 60),
#'                        female = c(0, 0, 1))
#' states <- data.frame(state_id = c(1, 2))
#' hesim_dat <- hesim_data(strategies = strategies,
#'                         patients = patients,
#'                         states = states)
#'
#' # Utility varies by health state and patient group
#' utility_tbl <- stateval_tbl(data.frame(state_id = rep(states$state_id, 2),
#'                                        grp_id = rep(rep(c(1, 2)), each = nrow(states)), 
#'                                        mean = c(.8, .7, .75, .55),
#'                                        se = c(.18, .12, .10, .06)),
#'                             dist = "beta",
#'                             hesim_data = hesim_dat)
#' print(utility_tbl)
#' utilmod <- create_StateVals(utility_tbl, n = 2)
#'
#' # Costs vary by treatment strategy
#' cost_tbl <- stateval_tbl(data.frame(strategy_id = strategies$strategy_id,
#'                                     mean = c(5000, 3000),
#'                                     se = c(200, 100)),
#'                          dist = "gamma",
#'                          hesim_data = hesim_dat)
#' print(cost_tbl)
#' costmod <- create_StateVals(cost_tbl, n = 2)
#'
#'
#' @export
stateval_tbl <- function(tbl, dist = c("norm", "beta", "gamma", 
                                       "lnorm", "unif", "fixed", "custom"),
                         hesim_data = NULL){
  x <- ParamsTbl$new(data.table(tbl), match.arg(dist), hesim_data)
  x$create_tbl("stateval_tbl")
  return(x$tbl[, ])
}

# StateVals --------------------------------------------------------------------
#' Create \code{StateVals} object
#' 
#' \code{create_StateVals} is a generic function for creating an object of class
#'  \code{\link{StateVals}} from a fitted statistical model or a \code{\link{stateval_tbl}}
#'  object. 
#' @param object A model object of the appropriate class.
#' @param input_data An object of class "expanded_hesim_data" returned by 
#' \code{\link{expand.hesim_data}}. Must be expanded by the data tables "strategies",
#' "patients", and "states".
#' @param n Number of random observations of the parameters to draw when parameters 
#' are fit using a statistical model.
#' @param point_estimate If \code{TRUE}, then the point estimates are returned and and no samples are drawn.
#' @param time_reset If \code{TRUE}, then time intervals reset each time a patient enters a new health 
#' state. See \code{\link{input_mats}}.
#' @param ... Further arguments passed to or from other methods. Currently unused. 
#' @return Returns an \code{\link{R6Class}} object of class \code{\link{StateVals}}.
#' @seealso \code{\link{StateVals}}
#' @export
create_StateVals <- function(object, ...){
  UseMethod("create_StateVals", object)
} 
 
#' @rdname create_StateVals
#' @export  
create_StateVals.lm <- function(object, input_data = NULL, n = 1000,
                                point_estimate = FALSE, ...){
  params <- create_params(object, n, point_estimate) 
  input_mats <- create_input_mats(object, input_data)
  return(StateVals$new(params = params, input_mats = input_mats))
}

#' @rdname create_StateVals 
#' @export
create_StateVals.stateval_tbl <- function(object, n = 1000, time_reset = FALSE, ...){
  x <- CreateFromParamsTbl$new(object, n)
  x$prep()
  stateval <- StateVals$new(params = x$params)
  stateval$params$time_reset <- time_reset
  return(stateval)
}


# Manual documentation in StateVals.Rd
#' @export
StateVals <- R6::R6Class("StateVals",
  public = list(
    input_mats = NULL,
    params = NULL,

    initialize = function(params, input_mats = NULL, time_reset = FALSE) {
      self$params <- params
      self$input_mats <- input_mats
    },
    
    sim = function(t, type = c("predict", "random")){
      type <- match.arg(type)
      self$check()
      res <- data.table(C_statevals_sim(self, t, type))
      res[, sample := sample + 1]
      return(res[])
    },
    
    check = function(){
      if(!inherits(self$params, c("tparams_mean", "params_lm"))){
        stop("Class of 'params' is not supported. See documentation.",
             call. = FALSE)
      }      
      if(!inherits(self$input_mats, c("input_mats", "NULL"))){
        stop("'input_mats' must be an object of class 'input_mats'",
            call. = FALSE)
      }
    }
  )
)

