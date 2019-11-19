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
  dist <- match.arg(dist)
  tbl <- data.table(tbl)
  tbl2 <- copy(tbl)
  cols <- colnames(tbl2)
  
  # Time intervals in the correct format
  if (!is.null(tbl2$time_start)){
    time_intervals <- time_intervals(unique(tbl2$time_start)) 
    pos <- match(tbl2$time_start, time_intervals$time_start)
    tbl2[, "time_id" := time_intervals$time_id[pos]]
    tbl2[, "time_stop" := time_intervals$time_stop[pos]]
  }
  
  # Check
  ## Column names
  check_column <- function(var){
    if (is.null(tbl2[[var]])){
      name <- switch(var,
                     "state_id" = "states",
                     "strategy_id" = "strategies",
                     "grp_id" = "patients")
      if (is.null(hesim_data[[name]])){
        msg <- paste0("If '", var, "' is not a column in 'tbl' ",
                      "then 'hesim_data' must be included as an argument ",
                      "and '",  name, "' must be an element of 'hesim_data'.")
        stop(msg, call. = FALSE)
      }
    }
  }
  check_column("state_id")
  check_column("strategy_id")
  check_column("grp_id")
  
  ## Need sample
  if (dist != "custom") {
    if ("sample" %in% colnames(tbl2)){
      stop(paste0("If 'sample' is in 'tbl', then 'dist' must equal 'custom'."))
    }
  }
  
  ## Distributions 
  if (dist == "norm"){
    if (!all(c("mean", "sd") %in% cols)){
      msg <- stop("If a normal distribution is specified, then tbl must ",
                  "contain the columns 'mean' and 'sd'.")
      stop(msg, call. = FALSE)         
    }
  } else if (dist == "beta"){
    if (!all(c("mean", "se") %in% cols) &
        !all(c("shape1", "shape2") %in% cols)){
      msg <- stop("If a beta distribution is specified, then tbl must either ",
                  "contain the columns 'mean' and 'se' or 'shape1' and 'shape2'.")
      stop(msg, call. = FALSE)      
    }
  } else if (dist == "gamma"){
    if (!all(c("mean", "se") %in% cols) &
        !all(c("shape", "rate") %in% cols) &
        !all(c("shape", "scale") %in% cols)){
      msg <- stop("If a gamma distribution is specified, then tbl must either ",
                  "contain the columns 'mean' and 'se', 'shape' and 'rate', ",
                  "or 'shape' and 'scale'.")
      stop(msg, call. = FALSE)        
    }
  } else if (dist == "lnorm"){
    if (!all(c("meanlog", "sdlog") %in% cols)){
      msg <- stop("If a lognormal distribution is specified, then tbl must ",
                  "contain the columns 'meanlog' and 'sdlog'.")
      stop(msg, call. = FALSE) 
    }
  } else if (dist == "unif"){
    if (!all(c("min", "max") %in% cols)){
      msg <- stop("If a uniform distribution is specified, then tbl must ",
                  "contain the columns 'min' and 'max'.")
      stop(msg, call. = FALSE)
    }
  } else if (dist == "fixed"){
    if (!all(c("est") %in% cols)){
      msg <- stop("If 'dist' = 'fixed', then tbl must ",
                  "contain the column 'est'.")
      stop(msg, call. = FALSE)
    }    
  } else if (dist == "custom"){
    if (!all(c("sample", "value") %in% cols)){
      msg <- stop("If 'dist' = 'custom', then tbl must ",
                  "contain the columns 'sample' and 'value'.")
      stop(msg, call. = FALSE)
    }  
  }
  
  ## Unique rows
  id_vars_all <- c("sample", "strategy_id", "state_id", "grp_id", "time_start")  
  id_vars <- id_vars_all[which(id_vars_all %in% colnames(tbl2))]
  if (length(id_vars) == 1){
    id_vars_msg <- id_vars
  } else{
    id_vars_msg <- paste0(paste(id_vars[1:length(id_vars) - 1], collapse = ", "),
                          ", and ", id_vars[length(id_vars)])
  }
  if (!all(tbl2[, .N, by = id_vars]$N == 1)) {
    stop(paste0("There must only be one row for each combination of ",
                id_vars_msg,
                " in 'tbl'."),
         call. = FALSE)
  }
  
  ## Number of rows
  size <- function(var){
    if (is.null(tbl2[[var]])){
      n <- 1
    } else{
      n <- length(unique(tbl2[[var]]))
    }
    return(n)
  }
  expected_n_samples <- size("sample")
  expected_n_strategies <- size("strategy_id")
  expected_n_states <- size("state_id")
  expected_n_grps <- size("grp_id")
  expected_n_times <- size("time_start")
  expected_n <- expected_n_samples * expected_n_strategies * expected_n_states * expected_n_grps * expected_n_times
  if (nrow(tbl2) != expected_n) {
    if (length(id_vars) == 1){
      stop(paste0("The number of rows in 'tbl' should equal ", expected_n, 
                  " which is the number of unique values of ",
                  id_vars_msg, " in 'tbl'")) 
    } else{
      stop(paste0("The number of rows in 'tbl' should equal ", expected_n, 
                  " which is the product of the number of unique values of ",
                  id_vars_msg, " in 'tbl'")) 
    }
  }
  
  # Nice sorting
  id_cols <- c("sample", "strategy_id", "patient_id", "grp_id", "state_id", "time_id",
               "time_start", "time_stop") 
  pos <- which(id_cols %in% colnames(tbl2))
  setcolorder(tbl2, id_cols[pos])
  
  # Return
  setattr(tbl2, "class", c("stateval_tbl", "data.table", "data.frame"))
  setattr(tbl2, "dist", dist)
  setattr(tbl2, "strategy_id", hesim_data$strategies$strategy_id)
  setattr(tbl2, "patients", data.table(hesim_data$patients))
  setattr(tbl2, "state_id", hesim_data$states$state_id)
  return(tbl2)
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
  
  # Random number generation
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
  
  # Transform
  if (attr(object, "dist") != "custom"){
    mu <- matrix(mu, ncol = n, byrow = FALSE) 
  } else {
    setorderv(tbl, "sample") 
    if (!is.null(n)){
      n_samples <- length(unique(tbl$sample))
      mu <- matrix(mu, ncol = n_samples, byrow = FALSE)
      samples <- sample_from_posterior(n = n, n_samples = n_samples)
      if (n != n_samples){
        mu <- mu[, samples, drop = FALSE]
      }
    }
    tbl <- tbl[sample == 1] 
  }
  tbl[, ("row_num") := 1:.N] 
  
  # Expand
  ## By strategy_id, state_id, and/or time interval
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
  
  ## By patient
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
  
  # Create object
  if (!is.null(tbl$time_id)){
    time_intervals <- unique(object[, c("time_id", "time_start", "time_stop")]) 
  } else{
    time_intervals <- NULL
  }
  
  tparams <- new_tparams_mean(value = mu,
                              n_samples = n,
                              strategy_id = tbl$strategy_id,
                              n_strategies = length(unique(tbl$strategy_id)),
                              patient_id = tbl$patient_id,
                              n_patients = length(unique(tbl$patient_id)),
                              state_id = tbl$state_id,
                              n_states = length(unique(tbl$state_id)),
                              time_id = tbl$time_id,
                              time_intervals = time_intervals,
                              n_times = nrow(time_intervals)) 
  return(StateVals$new(params = tparams, time_reset = time_reset))
}

create_StateVals.eval_model <- function(object, cost = TRUE, name = NULL,
                                        time_reset = FALSE, ...){
  out <- if (cost) object[["costs"]][[name]] else object$utility
  n_states <- object$n_states - 1 # The non-death states
  id  <- object$id[[attr(out, "id_index")]]
  out_id <- id[rep(1:nrow(id), each = n_states)]
  if ((is.numeric(out) && length(dim(out)) <= 1) || ncol(out) == 1){
    out_dt <- cbind(out_id, value = rep(out, each = n_states))
  } else{
    out_dt <- cbind(out_id, value = as.vector(t(out)))
  }
  out_dt[, ("state_id") := rep(1:n_states, times = nrow(id))]
  setnames(out_dt, "patient_id", "grp_id")
  return(create_StateVals(stateval_tbl(out_dt, dist = "custom"),
                          n = object$n))
}

# Manual documentation in StateVals.Rd
#' @export
StateVals <- R6::R6Class("StateVals",
  public = list(
    params = NULL,
    input_mats = NULL,
    time_reset = NULL,

    initialize = function(params, input_mats = NULL, time_reset = FALSE) {
      self$params <- params
      self$input_mats <- input_mats
      self$time_reset = time_reset
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
      stopifnot(is.logical(self$time_reset))
    }
  )
)

# Weighted length of stay ------------------------------------------------------
sim_wlos <- function (object, ...) {
  UseMethod("sim_wlos", object)
}

sim_wlos.stateprobs <- function(object, statevalmods, categories, dr = .03,
                                method = c("trapz", "riemann_left", "riemann_right")){
  method <- match.arg(method)
  state_id <- NULL
  
  # Checks
  ## State probabilities
  if(is.null(object)){
    stop("You must first simulate health state probabilities.",
         call. = FALSE)
  }
  
  # Discount rate
  check_dr(dr)
  
  # The state value models
  expected_samples <- max(object$sample)
  for (i in 1:length(statevalmods)){
    ## Number of samples
    if (statevalmods[[i]]$params$n_samples != expected_samples){
      msg <- paste0("Number of samples in each state value model must equal to ",
                    " the number of samples in the 'stateprobs' object, which is ",
                    expected_samples)
      stop(msg, call. = FALSE)
    }
 
    ## Number of states
    if(length(unique(object$state_id)) != get_id_object(statevalmods[[i]])$n_states + 1){
      msg <- paste0("The number of states in each 'StateVals' model ", 
                    "must be one less (since state values cannot be applied to the ",
                    "death state) than the number of states in 'stateprobs'.")
      stop(msg, call. = FALSE)
    }
  }
  
  # Simulate
  res <- data.table(C_sim_wlos(object[state_id != max(state_id)],
                               statevalmods, dr, categories,
                               unique(object$t),
                               method))
  res[, sample := sample + 1]
  return(res[])
} 

#' Weighted length of stay
#' 
#' Simulate weighted length of stay in order to compute costs and quality-adjusted
#' life-years (QALYs). 
#' 
#' @param object Objects of class \code{\link{stateprobs}} or \code{disprog}. 
#' @param utility_model A single object of class \code{\link{StateVals}} used
#' to simulate utility.
#' @param cost_models A list of objects of class \code{\link{StateVals}} used
#' to simulate costs.
#' @param dr Discount rate. 
#' @param method Method used to integrate state values when computing 
#' weighted length of stay. Options are \code{trapz} for the trapezoid rule,
#' \code{riemann_left} left for a left Riemann sum, and  
#' \code{riemann_right} right for a right Riemann sum.
#' @param lys If \code{TRUE}, then life-years are simulated in addition to 
#' QALYs. 
#' @return \code{sim_costs} and \code{sim_qalys} return objects of class
#' \code{costs} and \code{qalys}, respectively. 
#' @details 
#' Discounted costs and QALYs are calculated by integrating the "weighted" probability of being in each state. 
#' Weights are a function of the discount factor and the state value predicted using either the cost or QALY model. 
#' Mathematically, discounted costs and QALYs in health state \eqn{s} are calculated as,
#'
#'\deqn{\int_0^T w_h e^{-rt} P_h(t) dt },
#'
#' where for health state \eqn{h} and time {t}, \eqn{w_h} is the predicted cost 
#' or QALY weight, \eqn{r} is the discount rate, and \eqn{P_h(t)} is the 
#' probability of being in a given health state. The integral is calculated
#'  numerically from the points in \code{t_} using the approach selected 
#'  using the argument \code{method}.
#'
#' @export
#' @name sim_wlos
sim_qalys <- function(object, utility_model, dr, method, lys){
  utility_model$check()
  qalys <- sim_wlos(object,
                    list(utility_model),
                    "qalys",
                    dr,
                    method)
  setnames(qalys, "value", "qalys")
  setattr(qalys, "class", 
          c("qalys", "data.table", "data.frame"))
  return(qalys)
}

#' @export
#' @rdname sim_wlos
sim_costs <- function(object, cost_models, dr, method){
  if(!is.list(cost_models)){
    stop("'cost_models' must be a list", call. = FALSE)
  }
  for (i in 1:length(cost_models)){
    cost_models[[i]]$check()
  }
  if (is.null(names(cost_models))){
    categories <- paste0("Type ", seq(1, length(cost_models)))
  } else{
    categories <- names(cost_models)
  }   
  costs <- sim_wlos(object,
                     cost_models,
                     categories,
                     dr,
                     method)
  setnames(costs, "value", "costs")
  setattr(costs, "class", 
          c("costs", "data.table", "data.frame"))
  return(costs)
}

#' Costs object
#'
#' An object of class \code{costs} returned from methods 
#' \code{$sim_costs()} in model classes that stores simulated costs. 
#' 
#' @section Components:
#' A \code{costs} object inherits from \code{data.table} and contains
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
#' An object of class \code{qalys} returned from methods 
#' \code{$sim_qalys()} in model classes that stores simulated 
#' quality-adjusted life-years (QALYs).
#' 
#' @section Components:
#' A \code{qalys} object inherits from \code{data.table} and contains
#' the following columns:
#' 
#' \describe{
#'   \item{sample}{A random sample from the PSA.}
#'   \item{strategy_id}{The treatment strategy ID.}
#'   \item{patient_id}{The patient ID.}
#'   \item{state_id}{The health state ID.}
#'   \item{dr}{The rate used to discount QALYs.}
#'   \item{category}{A single category always equal to "qalys".}
#'   \item{costs}{The simulated values of QALYs.}
#' }
#'
#' @name qalys
NULL

# State probability object -----------------------------------------------------
#' State probability object
#'
#' An object of class \code{stateprobs} returned from methods 
#' \code{$sim_stateprobs()} in model classes. 
#' 
#' @section Components:
#' A \code{\link{stateprobs}} object inherits from \code{data.table} and contains
#' the following columns:
#' 
#' \describe{
#'   \item{sample}{A random sample from the PSA.}
#'   \item{strategy_id}{The treatment strategy ID.}
#'   \item{patient_id}{The patient ID.}
#'   \item{state_id}{The health state ID.}
#'   \item{t}{The time at which a state probability is computed.}
#'   \item{prob}{The probability of being in a given health state.}
#' }
#' 
#' When simulating individual-level models, the \code{patient_id} column is
#' not included as state probabilities are computed by averaging across patients.
#'
#'    
#' @name stateprobs
NULL