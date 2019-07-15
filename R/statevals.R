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
#' @seealso \code{\link{create_StateVals}}
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
  return(StateVals$new(input_mats = input_mats, params = params))
}

#' @rdname create_StateVals 
#' @export
create_StateVals.stateval_tbl <- function(object, n = 1000, time_reset = FALSE, ...){

  # Parameters
  tbl <- copy(object)
  n_rows <- nrow(tbl)
  if (attr(object, "dist") == "norm"){
    mu <- stats::rnorm(n * n_rows, mean = tbl$mean, sd = tbl$sd)
  } else if (attr(object, "dist") == "beta"){
    if (all(c("shape1", "shape2") %in% colnames(tbl))){
      mu <- stats::rbeta(n * n_rows, shape1 = tbl$shape1, shape2 = tbl$shape2)
    } else if (all(c("mean", "se") %in% colnames(tbl))){
      mom_params <- mom_beta(tbl$mean, tbl$se)
      mu <- stats::rbeta(n * n_rows, shape1 = mom_params$shape1, shape2 = mom_params$shape2) 
    } 
  } else if (attr(object, "dist") == "gamma"){
      if (all(c("shape", "rate") %in% colnames(tbl))){
        mu <- stats::rgamma(n * n_rows, shape = tbl$shape, rate = tbl$rate)
      } else if (all(c("shape", "scale") %in% colnames(tbl))){
        mu <- stats::rgamma(n * n_rows, shape = tbl$shape, scale = tbl$scale)
      } else if (all(c("mean", "se") %in% colnames(tbl))){
        mom_params <- mom_gamma(tbl$mean, tbl$se)
        mu <- stats::rgamma(n * n_rows, shape = mom_params$shape, scale = mom_params$scale) 
      } 
  } else if (attr(object, "dist") == "lnorm"){
      mu <- stats::rlnorm(n * n_rows, meanlog = tbl$meanlog, sdlog = tbl$sdlog)
  } else if (attr(object, "dist") == "unif"){
      mu <- stats::runif(n * n_rows, min = tbl$min, max = tbl$max) 
  } else if (attr(object, "dist") == "fixed"){
      mu <- rep(tbl$est, times = n)
  } else if (attr(object, "dist") == "custom"){
      mu <- tbl$value
  }
  
  if (attr(object, "dist") != "custom"){
    mu <- matrix(mu, ncol = n, byrow = FALSE) 
  } else {
    setorderv(tbl, "sample") 
    if (!is.null(n)){
      n_samples <- length(unique(tbl$sample))
      mu <- matrix(mu, ncol = n_samples, byrow = FALSE)
      if (n < n_samples){
        samples <- sample.int(n_samples, n, replace = FALSE) 
      } else if (n > n_samples) {
        warning("'n' is larger than the number of unique values of 'sample' in 'tbl'.")
        samples <- sample.int(n_samples, n, replace = TRUE) 
      }
      if (n != n_samples){
        mu <- mu[, samples, drop = FALSE]
      }
    }
    tbl <- tbl[sample == 1] 
  }
  tbl[, ("row_num") := 1:.N] 
  
  ## Expand by strategy_id, state_id, and/or time interval
  tbl_list <- list()
  id_vars <- c("strategy_id", "state_id", "time_start")
  i <- 1
  for (var in id_vars){
    if (is.null(tbl[[var]])){
      if (!is.null(attr(object, var))){
        tbl_i <- data.frame(tmp_var = attr(object, var))
        setnames(tbl_i, "tmp_var", var)
        tbl_list[[i]] <- tbl_i
        i <- i + 1
      }
    }
  }
  tbl_list <- c(list(data.frame(tbl)), tbl_list)
  tbl <- Reduce(function(...) merge(..., by = NULL), tbl_list)
  tbl <- data.table(tbl)
  
 ## Expand by patient
  merge <- TRUE
  if (is.null(tbl$grp_id)){ # If group ID is not specified
    tbl[, ("grp_id") := 1]
    patient_lookup <- data.table(patient_id = attr(object, "patients")$patient_id, 
                                 grp_id = 1)
  } else { # Else if group ID is specified
    if (is.null(attr(object, "patients")$grp_id)) { # If the patient lookup table does not exist
      setnames(tbl, "grp_id", "patient_id")
      merge <- FALSE
    } else{
        patient_lookup <- attr(object, "patients")[, c("patient_id", "grp_id"), 
                                                    with = FALSE] 
    }
  }
  if (merge){
    tbl <- merge(tbl, patient_lookup, by = c("grp_id"), allow.cartesian = TRUE,
                 sort = FALSE) 
  }
  if (is.null(tbl[["time_start"]])){
    setorderv(tbl, cols = c("strategy_id", "patient_id", "state_id")) 
  } else{
    setorderv(tbl, cols = c("strategy_id", "patient_id", "state_id", 
                            "time_id")) 
  }
  mu <- mu[tbl$row_num, , drop = FALSE]

  ## Create object
  params <- new_params_mean(mu = mu, 
                            sigma = rep(0, n),
                            n_samples = n)
  # Input matrices
  if (!is.null(tbl$time_id)){
    time_intervals <- unique(object[, c("time_id", "time_start", "time_stop")]) 
  } else{
    time_intervals <- NULL
  }
  
  input_mats <- new_input_mats(X = NULL,
                              strategy_id = tbl$strategy_id,
                              n_strategies = length(unique(tbl$strategy_id)),
                              patient_id = tbl$patient_id,
                              n_patients = length(unique(tbl$patient_id)),
                              state_id = tbl$state_id,
                              n_states = length(unique(tbl$state_id)),
                              time_id = tbl$time_id,
                              time_intervals = time_intervals,
                              time_reset = time_reset,
                              n_times = nrow(time_intervals)) 
  return(StateVals$new(input_mats = input_mats, params = params))
}


# Manual documentation in StateVals.Rd
#' @export
StateVals <- R6::R6Class("StateVals",
  public = list(
    input_mats = NULL,
    params = NULL,

    initialize = function(input_mats, params, time_reset = FALSE) {
      self$input_mats <- input_mats
      self$params <- params
    },
    
    sim = function(t, type = c("predict", "random")){
      type <- match.arg(type)
      self$check()
      res <- data.table(C_statevals_sim(self, t, type))
      res[, sample := sample + 1]
      return(res[])
    },
    
    check = function(){
      if(!inherits(self$input_mats, "input_mats")){
        stop("'input_mats' must be an object of class 'input_mats'",
            call. = FALSE)
      }
      if(!inherits(self$params, c("params_mean", "params_lm"))){
          stop("Class of 'params' is not supported. See documentation.",
               call. = FALSE)
      }
    }
  )
)

