# StateVals class --------------------------------------------------------------
#' Model for state values
#' 
#' @description
#' Simulate values (i.e., utility or costs) associated with health states in a 
#' state transition or partitioned survival model. 
#' 
#' @example man-roxygen/example-StateVals.R
#' @name StateVals
NULL

#' @rdname StateVals
#' @export
StateVals <- R6::R6Class(
  "StateVals",
  
  public = list(
    #' @field params Parameters for simulating state values. Currently supports
    #'  objects of class [`tparams_mean`] or [`params_lm`].   
    params = NULL,
    
    #' @field input_data  An object of class [input_mats]. Only used for
    #' [`params_lm`] objects.
    input_data = NULL,
    
    #' @field method The method used to simulate costs and 
    #' quality-adjusted life-years (QALYs) as a function of state values.
    #'  If `wlos`, then costs and QALYs are
    #' simulated by weighting state values by the length of stay in a health
    #' state. If `starting`, then state values represent a one-time value
    #' that occurs when a patient enters a health state. When `starting` is 
    #' used in a cohort model, the state values only accrue at time 0; 
    #' in contrast, in an individual-level model, state values
    #' accrue each time a patient enters a new state and are discounted based on
    #' time of entrance into that state. When `transition`, then state values
    #' represent a value that occurs at the transition time; only used for
    #' individual-level models.
    method = NULL,
    
    #' @field time_reset If `FALSE` then time intervals are based on time since
    #'  the start of the simulation. If `TRUE`, then time intervals reset each 
    #'  time a patient enters a new health state. This is relevant if, for example, 
    #'  costs vary over time within health states. Only used if `method = wlos`.
    time_reset = NULL,    
    
    #' @description
    #' Create a new `StateVals` object.
    #' @param params The `params` field.
    #' @param input_data The `input_data` field.
    #' @param method The `method` field.
    #' @param time_reset The `time_reset` field.
    #' @return A new `StateVals` object.
    initialize = function(params, input_data = NULL,
                          method = c("wlos", "starting", "transition"),
                          time_reset = FALSE) {
      self$params <- params
      self$input_data <- input_data
      self$method <- match.arg(method)
      self$time_reset <- time_reset
    },
    
    #' @description
    #' Simulate state values with either predicted means or random samples by
    #'  treatment strategy, patient, health state, and time `t`.
    #' @param t A numeric vector of times. 
    #' @param type  `"predict"` for mean values or `"random"` for random samples. 
    #' @return A `data.table` of simulated state values with columns for `sample`,
    #' `strategy_id`, `patient_id`, `state_id`, `time`, and `value`.  
    sim = function(t, type = c("predict", "random")){
      type <- match.arg(type)
      self$check()
      res <- data.table(C_statevals_sim(self, sort(t), type))
      res[, sample := sample + 1]
      return(res[])
    },
    
    #' @description
    #' Input validation for class. Checks that fields are the correct type. 
    check = function(){
      if(!inherits(self$params, c("tparams_mean", "params_lm"))){
        stop("Class of 'params' is not supported. See documentation.",
             call. = FALSE)
      }      
      if(!inherits(self$input_data, c("input_mats", "NULL"))){
        stop("'input_data' must be an object of class 'input_mats'",
             call. = FALSE)
      }
      stopifnot(is.logical(self$time_reset))
    }
  )
)

# stateval_tbl -----------------------------------------------------------------
#' Table to store state value parameters
#' 
#' Create a table for storing parameter estimates used to simulate costs or 
#' utility in an economic model by treatment strategy, patient, health state, and
#' (optionally) time interval. 
#' 
#' @param tbl A `data.frame` or `data.table` for storing parameter 
#' values. See "Details" for specifics. 
#' @param dist Probability distribution used to sample parameters for a 
#' probabilistic sensitivity analysis (PSA). 
#' @param hesim_data A [`hesim_data`] object. This argument is deprecated
#' and should be passed to [`create_StateVals.stateval_tbl()`] instead.
#' @param grp_var The name of the variable used to group patients.
#' 
#' @details 
#' `tbl` is a tabular object containing columns for treatment 
#' strategies (`strategy_id`), patients (`patient_id`),
#' health states (`state_id`), and/or the start of time intervals 
#' (`time_start`). The table must contain at least one column
#' named `strategy_id`, `patient_id`, or `state_id`, 
#' but does not need to contain all of them. Each row denotes a unique 
#' treatment strategy, patient, health state, and/or time interval pair.
#' `tbl` may also contain a column with the name specified in `grp_var` 
#' (rather than `patient_id`) so that state values are assigned to 
#' groups of patients.
#' 
#' `tbl` must also contain columns summarizing the state values for each
#' row, which depend on the probability distribution selected with `dist`. 
#' Available distributions include the normal (`norm`), beta (`beta`),
#' gamma (`gamma`), lognormal (`lnorm`), and uniform (`unif`)
#'  distributions. In addition, the option `fixed` can be used if estimates
#'  are known with certainty and `custom` can be used if 
#'  parameter values for a PSA  have been previously
#' sampled from an arbitrary probability distribution.
#'  The columns in `tbl` that must be included,
#'  by distribution, are:
#' 
#' \describe{
#' \item{norm}{`mean` and `sd`}
#' \item{beta}{`mean` and `se` or `shape1` and `shape2`}
#' \item{gamma}{`mean` and `se`, `shape` and `rate`, 
#' or `shape` and `scale`}
#' \item{lnorm}{`meanlog` or `sdlog`}
#' \item{unif}{`min` and `max`}
#' \item{fixed}{`est`}
#' \item{custom}{`sample` and `value`}
#' }
#' 
#' Note that if `dist = "custom"`, then `tbl` must include a column 
#' named `sample` (an integer vector denoting a unique random draw) and
#'  `value` (denoting the value of the randomly sampled parameter). In this case, 
#'  there is a unique row in `tbl` for each random draw (`sample`) and
#'  each combination of strategies, patients, health states, and/or time intervals.
#' Again, `tbl` must contain at least one column
#' named `strategy_id`, `patient_id` (or `grp_var`), or `state_id`,
#'  but does not need to contain them all.
#'  
#'  
#' @return An object of class `stateval_tbl`, which is a `data.table` of
#' parameter values with attributes for `dist` and `grp_var`.
#' @seealso [`create_StateVals()`], [`StateVals`]
#' @examples 
#' strategies <- data.frame(strategy_id = c(1, 2))
#' patients <- data.frame(patient_id = seq(1, 3),
#'                        grp = c(1, 1, 2),
#'                        age = c(45, 50, 60),
#'                        female = c(0, 0, 1))
#' states <- data.frame(state_id = c(1, 2))
#' hesim_dat <- hesim_data(strategies = strategies,
#'                         patients = patients,
#'                         states = states)
#'
#' # Utility varies by health state and patient group
#' utility_tbl <- stateval_tbl(data.frame(state_id = rep(states$state_id, 2),
#'                                        grp = rep(rep(c(1, 2)), each = nrow(states)), 
#'                                        mean = c(.8, .7, .75, .55),
#'                                        se = c(.18, .12, .10, .06)),
#'                             dist = "beta",
#'                             grp_var = "grp")
#' print(utility_tbl)
#' utilmod <- create_StateVals(utility_tbl, n = 2, hesim_data = hesim_dat)
#'
#' # Costs vary by treatment strategy
#' cost_tbl <- stateval_tbl(data.frame(strategy_id = strategies$strategy_id,
#'                                     mean = c(5000, 3000),
#'                                     se = c(200, 100)),
#'                          dist = "gamma")
#' print(cost_tbl)
#' costmod <- create_StateVals(cost_tbl, n = 2, hesim_data = hesim_dat)
#'
#'
#' @export
stateval_tbl <- function(tbl, dist = c("norm", "beta", "gamma", 
                                       "lnorm", "unif", "fixed", "custom"),
                         hesim_data = NULL,
                         grp_var = NULL){
  if (!missing("hesim_data")) {
    warning("'hesim_data' argument deprecated; pass to create_StateVals() instead.")
  }
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
  id_vars_all <- c("sample", "strategy_id", "state_id", "patient_id", grp_var, "time_start")  
  id_vars <- id_vars_all[which(id_vars_all %in% colnames(tbl2))]
  if (length(id_vars) == 1){
    id_vars_msg <- id_vars
  } else if (length(id_vars) == 2){
    id_vars_msg <- paste0(id_vars[1], " and ", id_vars[2])
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
    if (is.null(var) || is.null(tbl2[[var]])){
      n <- 1
    } else{
      n <- length(unique(tbl2[[var]]))
    }
    return(n)
  }
  expected_n_samples <- size("sample")
  expected_n_strategies <- size("strategy_id")
  expected_n_states <- size("state_id")
  expected_n_patients <- size("patient_id")
  if (!is.null(tbl2[["patient_id"]])){
    expected_n_grps <- expected_n_patients
  } else{
    expected_n_grps <- size(grp_var)
  }
  expected_n_times <- size("time_start")
  expected_n <- expected_n_samples * expected_n_strategies * expected_n_states * expected_n_grps * expected_n_times
  if (nrow(tbl2) != expected_n) {
    stop(paste0("The number of rows in 'object' should equal ", expected_n, 
                " which is the number of unique values of ",
                id_vars_msg, " in 'object'."))
  }
  
  # Nice sorting
  id_cols <- c("sample", "strategy_id", "patient_id", grp_var, "state_id", "time_id",
               "time_start", "time_stop") 
  pos <- which(id_cols %in% colnames(tbl2))
  setcolorder(tbl2, id_cols[pos])
  
  # Return
  setattr(tbl2, "class", c("stateval_tbl", "data.table", "data.frame"))
  setattr(tbl2, "dist", dist)
  setattr(tbl2, "grp_var", grp_var)
  if (!is.null(hesim_data)) {
    attr(tbl2, "hesim_data") <- hesim_data
  }
  return(tbl2)
}

# create_StateVals methods -----------------------------------------------------
#' Create a `StateVals` object
#' 
#' `create_StateVals()` is a generic function for creating an object of class
#'  [`StateVals`] from a fitted statistical model or a [`stateval_tbl`]
#'  object. 
#' @param object A model object of the appropriate class.
#' @param input_data An object of class [`expanded_hesim_data`][expand.hesim_data()].
#' Must be expanded by treatment strategies, patients, and health states.
#' @param hesim_data A [`hesim_data`] object. Only required when `object` is of class
#' [`stateval_tbl`]. See "details". 
#' @param n Number of random observations of the parameters to draw when parameters 
#' are fit using a statistical model.
#' @param uncertainty Method determining how parameter uncertainty should be handled. See
#'  documentation in [`create_params()`].
#' @param ... Further arguments (`time_reset` and `method`) passed to [`StateVals$new()`][StateVals].
#' @details If `object` is a `stateval_tbl`, then a [`hesim_data`] object is used
#'  to specify treatment strategies, patients, and/or health states not included as 
#'  columns in the table, or, to match patients in the table to groups. Not required if 
#'  the table includes one row for each treatment strategy, patient, and health state
#'   combination. Patients are matched to groups by specifying both a `patient_id` 
#'   and a `grp_var` column in the `patients` table.
#' @return A [`StateVals`] object.
#' @seealso See [`StateVals`] for documentation of the class and additional examples. 
#' An example use case for [create_StateVals.stateval_tbl()] is provided in 
#' the [stateval_tbl()] documentation.
#' @export
create_StateVals <- function(object, ...){
  UseMethod("create_StateVals", object)
} 
 
#' @rdname create_StateVals
#' @example man-roxygen/example-create_StateVals.lm.R
#' @export  
create_StateVals.lm <- function(object, input_data = NULL, n = 1000,
                                uncertainty = c("normal", "none"), ...){
  uncertainty <- deprecate_point_estimate(list(...)$point_estimate, uncertainty,
                                          missing(uncertainty))
  params <- create_params(object, n, uncertainty) 
  input_mats <- create_input_mats(object, input_data)
  return(StateVals$new(params = params, input_data = input_mats, ...))
}

#' @rdname create_StateVals
#' @param method String matching a list of methods used for the `StateVals` object.
#' @export
create_StateVals.stateval_tbl <- function(object, hesim_data = NULL, n = 1000,
                                          method = c("wlos","starting","transition"), ...){

  method <- match.arg(method)
  grp_var <- attr(object, "grp_var")
  
  # For backwards compatibility, use hesim_data attribute of object 
  # if needed and available
  if (!is.null(attr(object, "hesim_data")) && is.null(hesim_data)) {
    hesim_data <- attr(object, "hesim_data")
  } 
  
  # Check whether hesim_data is required
  need_hesim_data <- function(var){
    if (is.null(object[[var]])){
      name <- switch(var,
                     "state_id" = "states",
                     "strategy_id" = "strategies",
                     "patient_id" = "patients")
      if (is.null(hesim_data[[name]])){
        msg <- paste0("If '", var, "' is not a column in 'object', ",
                      "then 'hesim_data' must be included as an argument ",
                      "and '",  name, "' must be an element of 'hesim_data'.")
        stop(msg, call. = FALSE)
      }
    }
  }
  need_hesim_data("state_id")
  need_hesim_data("strategy_id")
  need_hesim_data("patient_id")
  
  ## Additional logic for patients
  if (is.null(object[["patient_id"]])){
    if (is.null(hesim_data[["patients"]][["patient_id"]]) & 
        (is.null(grp_var) || is.null(object[[grp_var]]))){
      stop(paste0("If 'patient_id' is not included as a column in `object`, ",
                  "then both 'patient_id' and 'grp_var' cannot be missing from the ",
                  " 'patients' element of hesim_data"),
           call. = FALSE)
    }
  }

  # If not NULL, add hesim_data to attributes of object 
  if (!is.null(hesim_data)) {
    setattr(object, "strategy_id", hesim_data$strategies$strategy_id)
    setattr(object, "patients", data.table(hesim_data$patients))
    setattr(object, "state_id", hesim_data$states$state_id)
  }
  
  # Random number generation
  tbl <- copy(object)
  n_rows <- nrow(tbl)
  if (attr(object, "dist") == "norm"){
    mu <- stats::rnorm(n * n_rows, mean = tbl$mean, sd = tbl$sd)
  } else if (attr(object, "dist") == "beta"){
    if (all(c("shape1", "shape2") %in% colnames(tbl))){
      mu <- stats::rbeta(n * n_rows, shape1 = tbl$shape1, shape2 = tbl$shape2)
    } else if (all(c("mean", "se") %in% colnames(tbl))){
      mu <- mom_fun_rng(n, rng_fun = "rbeta", mom_fun = "mom_beta",
                        mean = tbl$mean, sd = tbl$se)
    } 
  } else if (attr(object, "dist") == "gamma"){
    if (all(c("shape", "rate") %in% colnames(tbl))){
      mu <- stats::rgamma(n * n_rows, shape = tbl$shape, rate = tbl$rate)
    } else if (all(c("shape", "scale") %in% colnames(tbl))){
      mu <- stats::rgamma(n * n_rows, shape = tbl$shape, scale = tbl$scale)
    } else if (all(c("mean", "se") %in% colnames(tbl))){
      mu <- mom_fun_rng(n, rng_fun = "rgamma", mom_fun = "mom_gamma",
                        mean = tbl$mean, sd = tbl$se)
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
  grp_var <- attr(object, "grp_var")
  patients <- copy(attr(object, "patients"))
  if (is.null(tbl[["patient_id"]])){
    if (is.null(grp_var)){
      grp_var <- "grp"
      tbl[, ("grp") := 1]
      patients$grp <- 1
    }
    tbl <- merge(tbl, patients, by = grp_var, 
                 allow.cartesian = TRUE, sort = FALSE) 
  }
  
  # Sorting
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
  transp <- (method == "transition")
  tparams <- new_tparams_mean(value = mu,
                              n_samples = n,
                              strategy_id = tbl$strategy_id,
                              n_strategies = length(unique(tbl$strategy_id)),
                              patient_id = tbl$patient_id,
                              n_patients = length(unique(tbl$patient_id)),
                              state_id = if (transp) NULL else tbl$state_id,
                              transition_id = if (transp) tbl$state_id else NULL,
                              n_states = if (transp) NULL else length(unique(tbl$state_id)),
                              n_transitions = if (transp) length(unique(tbl$state_id)) else NULL,
                              time_id = tbl$time_id,
                              time_intervals = time_intervals,
                              n_times = nrow(time_intervals),
                              grp_id = tbl$grp_id,
                              patient_wt = tbl$patient_wt) 
  return(StateVals$new(params = tparams, method=method, ...))
}

#' @export
create_StateVals.eval_model <- function(object, cost = TRUE, name = NULL,
                                        init_args = NULL, ...){
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
  return(do.call("create_StateVals", 
                 args = c(list(object = stateval_tbl(out_dt, dist = "custom"),
                               n = object$n),
                           init_args)))
}

