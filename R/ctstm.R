# CtstmTrans -------------------------------------------------------------------
#' An `R6` base class for continuous time state transition models
#'
#' @description
#' A non-exported base class for continuous time state transition models containing methods
#' that can be used to summarize fitted multi-state models. 
#' @format An [R6::R6Class] object.
#' @seealso [create_IndivCtstmTrans()], [IndivCtstmTrans]
CtstmTrans <- R6::R6Class("CtstmTrans",
  private = list(
    check_base = function(){
      if(!inherits(self$input_data, "input_mats")){
        stop("'input_data' must be an object of class 'input_mats'",
            call. = FALSE)
      }
      if(!inherits(self$params, c("params_surv", 
                                  "params_surv_list"))){
        stop("Class of 'params' is not supported. See documentation.",
             call. = FALSE)
      }
      if (!is.matrix(self$trans_mat)){
        stop("'trans_mat' must be a matrix",
             call. = FALSE)
      }
      if(nrow(self$trans_mat) != ncol(self$trans_mat)){
        stop("'trans_mat' must be a square matrix",
             call. = FALSE)
      }
      if(inherits(self$params, "params_surv_list")){
        if(max(self$trans_mat, na.rm = TRUE) != length(self$params)){
          stop(paste0("The number of models in 'params' must equal the maximum integer in ",
                      "in 'trans_mat'."),
               call. = FALSE)
        }
      }
    },
    
      summary = function(t, type = c("hazard", "cumhazard")){
        self$check()
        type <- match.arg(type)
        res <- data.table(C_ctstm_summary(self, t, type))
        res[, trans := trans + 1]
        res[, sample := sample + 1]
        if (type == "hazard") setnames(res, "value", "hazard")
        if (type == "cumhazard") setnames(res, "value", "cumhazard")
        return(res[])
      }
    ), # end private
  
  public = list(
    #' @description
    #' Predict the hazard functions for each health state transition.
    #' @param t  A numeric vector of times.
    #' @return A `data.table` with columns `trans` (the transition number),
    #'  `sample`, `strategy_id`, `grp_id`, `t`, and `hazard`.
    hazard = function(t){
      private$summary(t, "hazard")
    },

    #' @description
    #' Predict the cumulative hazard functions for each health state transition.
    #' @param t  A numeric vector of times.    
    #' @return A `data.table` with columns `trans`,
    #'  `sample`, `strategy_id`, `grp_id`, `t`, and `hazard`.
    cumhazard = function(t){
       private$summary(t, "cumhazard")
    }
  )
)

# IndivCtstmTrans --------------------------------------------------------------
indiv_ctstm_sim_disease <- function(trans_model, max_t = 100, max_age = 100,
                                    progress = NULL){
  sample <- from <- to <- NULL # to avoid no visible bindings CRAN warning
  if (any(trans_model$start_age > max_age)){
    stop("Starting ages in the simulation must be less than maximum age.",
         call. = FALSE)
  }
  if (is.null(progress)){
    progress <- 0
  }
  
  # Simulate
  disprog <- C_ctstm_sim_disease(trans_model, trans_model$start_state - 1, 
                                 trans_model$start_age,
                                 rep(0, trans_model$input_data$n_patients), # Start time is always 0
                                 trans_model$death_state - 1, 
                                 trans_model$clock, trans_model$reset_states - 1,
                                 max_t, max_age, progress)
  disprog <- data.table(disprog)
  disprog[, sample := sample + 1]
  disprog[, from := from + 1]
  disprog[, to := to + 1]
  return(disprog[, ])
}

indiv_ctstm_sim_stateprobs <- function(disprog = NULL, trans_model = NULL, t, ...){
  
  # Simulate disease progression if 'disprog' is missing
  if (is.null(disprog)){
    trans_model$check()
    dots <- list(...)
    fun <- trans_model$sim_disease
    disprog <- do.call("fun", dots)
    
  } else{
    disprog <- copy(disprog)
  }
  
  # Compute state probabilities
  ## Indexing for C++
  disprog[, strategy_index := .GRP, by = "strategy_id"]
  strategy_index <- disprog$strategy_index - 1
  disprog[, strategy_index := NULL]
  
  disprog[, grp_index := .GRP, by = "grp_id"]
  grp_index <- disprog$grp_index - 1
  disprog[, grp_index := NULL]

  ## Dimensions of simulation
  n_grps <- length(unique(disprog$grp_id))
  n_states <- nrow(trans_model$trans_mat)
  n_strategies <- trans_model$input_data$n_strategies
  n_patients <- trans_model$input_data$n_patients
  if (inherits(trans_model$params, "params_surv_list")){
    n_samples <- trans_model$params[[1]]$n_samples
  } else{
    n_samples <- trans_model$params$n_samples
  }
  unique_strategy_id <- unique(trans_model$input_data$strategy_id)
  unique_grp_id <- unique(trans_model$input_data$grp_id)
  
  ## Computation
  stprobs <- C_ctstm_indiv_stateprobs(disprog, 
                                      t, 
                                      n_samples,
                                      n_strategies, 
                                      unique_strategy_id,
                                      strategy_index,
                                      n_grps,
                                      unique_grp_id,
                                      grp_index,
                                      n_states,
                                      n_patients)
  stprobs <- data.table(stprobs)
      
  ## C++ to R indexing
  stprobs[, "sample" := get("sample") + 1]
  stprobs[, "state_id" := get("state_id") + 1]
  
  # Return
  return(stprobs[])      
}

#' Create \code{IndivCtstmTrans} object
#' 
#' A generic function for creating an object of class \code{\link{IndivCtstmTrans}}.
#' @param object A fitted survival model or the parameters of a survival model.  
#' @param input_data An object of class "expanded_hesim_data" returned by 
#' \code{\link{expand.hesim_data}}.
#' @param n Number of random observations of the parameters to draw.
#' @param trans_mat The transition matrix describing the states and transitions in a 
#' multi-state model in the format from the \link[mstate]{mstate} package. See \link{IndivCtstmTrans}.
#' @param clock "reset" for a clock-reset model, "forward" for a clock-forward model, and "mix" for a mixture
#' of clock-reset and clock-forward models. See the field \code{clock} in \code{\link{IndivCtstmTrans}}.
#' @param reset_states A vector denoting the states in which time resets. See the field 
#' \code{reset_states} in \code{\link{IndivCtstmTrans}}.
#' @param point_estimate If \code{TRUE}, then the point estimates are returned and and no samples are drawn.
#' @param ... Further arguments passed to \code{IndivCtstmTrans$new()} in \code{\link{IndivCtstmTrans}}.
#' @return Returns an \code{\link{R6Class}} object of class \code{\link{IndivCtstmTrans}}.
#' @seealso \code{\link{IndivCtstmTrans}}
#' @name create_IndivCtstmTrans
#' @rdname create_IndivCtstmTrans
#' @export
create_IndivCtstmTrans <- function(object, ...){
  UseMethod("create_IndivCtstmTrans", object)
}

#' @export
#' @rdname create_IndivCtstmTrans
create_IndivCtstmTrans.flexsurvreg_list <- function(object, input_data, trans_mat, clock = c("reset", "forward"),
                                                    n = 1000, point_estimate = FALSE, ...){
  input_mats <- create_input_mats(object, input_data)
  params <- create_params(object, n = n, point_estimate = point_estimate)
  return(IndivCtstmTrans$new(input_data = input_mats, params = params, trans_mat = trans_mat,
                             clock = match.arg(clock), ...))
}

#' @export
#' @rdname create_IndivCtstmTrans
create_IndivCtstmTrans.flexsurvreg <- function(object, input_data, trans_mat, clock = c("reset", "forward"),
                                               n = 1000, point_estimate = FALSE, ...){
  input_mats <- create_input_mats(object, input_data)
  params <- create_params(object, n = n, point_estimate = point_estimate)
  return(IndivCtstmTrans$new(input_data = input_mats, params = params, trans_mat = trans_mat, 
                             clock = match.arg(clock), ...))
}

#' @export
#' @rdname create_IndivCtstmTrans
create_IndivCtstmTrans.params_surv <- function(object, input_data, trans_mat, 
                                               clock = c("reset", "forward", "mix"),
                                               reset_states = NULL,...){
  input_mats <- create_input_mats(object, input_data)
  return(IndivCtstmTrans$new(input_data = input_mats, params = object, trans_mat = trans_mat,
                             clock = match.arg(clock), reset_states = reset_states, ...))
}

#' Transitions for an individual-level continuous time state transition model
#'
#' @description
#' Simulate health state transitions in an individual-level continuous time state
#'  transition model with parameters that were estimated using a multi-state model.
#' @format An [R6::R6Class] object.
#' @examples 
#' library("flexsurv")
#'
#' # Simulation data
#' strategies <- data.frame(strategy_id = c(1, 2, 3))
#' patients <- data.frame(patient_id = seq(1, 3),
#'                        age = c(45, 50, 60),
#'                        female = c(0, 0, 1))
#'
#' # Multi-state model with transition specific models
#' tmat <- rbind(c(NA, 1, 2),
#'               c(NA, NA, 3),
#'               c(NA, NA, NA))
#' fits <- vector(length = max(tmat, na.rm = TRUE), mode = "list")
#' for (i in 1:length(fits)){
#'   fits[[i]] <- flexsurvreg(Surv(years, status) ~ 1, 
#'                            data = bosms3[bosms3$trans == i, ],
#'                            dist = "exp")
#' }
#' fits <- flexsurvreg_list(fits)
#'
#' # Simulation model
#' hesim_dat <- hesim_data(strategies = strategies,
#'                         patients = patients)
#' fits_data <- expand(hesim_dat)
#' transmod <- create_IndivCtstmTrans(fits, input_data = fits_data, 
#'                                    trans_mat = tmat,
#'                                    n = 2,
#'                                    point_estimate = FALSE)       
#' head(transmod$hazard(c(1, 2, 3)))
#' head(transmod$cumhazard(c(1, 2, 3)))
#' transmod$sim_stateprobs(t = c(0, 5, 10))[t == 5]
#' @seealso [create_IndivCtstmTrans()], [IndivCtstm]
#' @export
IndivCtstmTrans <- R6::R6Class("IndivCtstmTrans",
  inherit = CtstmTrans,
  
  private = list(
    
    check_history = function(field){
      field_name <- deparse(substitute(field))
      if (length(field) !=1 & length(field) != self$input_data$n_patients){
        stop(paste0("The length of '", field_name, "' must either be 1 or the number ",
                          "of simulated patients."),
                 call. = FALSE)        
      }
      if (length(field) == 1){
        field <- rep(field, self$input_data$n_patients)
      }
      return(field)
    }
  ), # end private

  public = list(
    #' @field params An object of class [params_surv] or [params_surv_list].
    params = NULL,
    
    #' @field input_data Input data used to simulate health state transitions 
    #' by sample from the probabilistic sensitivity analysis (PSA), treatment strategy and patient.
    #' Must be an object of class [input_mats]. If `params` contains parameters from
    #' a list of models (i.e., of class [params_surv_list]), then `input_data` 
    #' must contain a unique row for each treatment strategy
    #' and patient; if `params` contains parameters from a joint model 
    #' (i.e., of class [params_surv]), then `input_data` must contain a unique
    #' row for each treatment strategy, patient, and transition.
    input_data = NULL,
    
    #' @field trans_mat A transition matrix describing the states and transitions 
    #' in a multi-state model in the format from the [mstate][mstate::mstate] package. 
    #' See the documentation for the argument `"trans"` in [mstate::msprep].
    trans_mat = NULL, 
    
    #' @field start_state A scalar or vector denoting the starting health state. 
    #' Default is the first health state. If a vector, must be equal to the number of simulated patients.
    start_state = NULL,
    
    #' @field start_age A scalar or vector denoting the starting age of each patient 
    #' in the simulation. Default is 38. If a vector, must be equal to the number of simulated patients.
    start_age = NULL,
    
    #' @field death_state The death state in `trans_mat`. Used with `max_age` 
    #' in `sim_disease` as patients transition to this state upon reaching maximum age.
    #'  By default, it is set to the final absorbing state (i.e., a row in `trans_mat` with all NAs).
    death_state = NULL,
    
    #' @field clock "reset" for a clock-reset model, "forward" for a clock-forward model, 
    #' and "mix" for a mixture of clock-reset and clock-forward models. A clock-reset model 
    #' is a semi-Markov model in which transition rates depend on time since entering a state. 
    #' A clock-forward model is a Markov model in which transition rates depend on time 
    #' since entering the initial state. If `"mix"` is used, then 
    #' `reset_states` must be specified.
    clock = NULL,
    
    #' @field reset_states A vector denoting the states in which time resets. 
    #' Hazard functions are always a function of elapsed time since either the 
    #' start of the model or from when time was previously reset. Only used if 
    #' `clock = "mix"`.
    reset_states = NULL,
    
    #' @description
    #' Create a new `IndivCtstmTrans` object.
    #' @param params The `params` field.
    #' @param input_data The `input_data` field.
    #' @param trans_mat The `trans_mat` field.
    #' @param start_state The `start_state` field.
    #' @param start_age The `start_age` field.
    #' @param death_state The `death_state` field.
    #' @param clock The `clock` field.
    #' @param reset_states The `reset_states` field.
    #' @return A new `IndivCtstmTrans` object.    
    initialize = function(params, input_data, trans_mat, 
                          start_state = 1,
                          start_age = 38,
                          death_state = NULL,
                          clock = c("reset", "forward", "mix"),
                          reset_states = NULL) {
      self$params <- params
      self$input_data <- input_data
      self$trans_mat <- trans_mat
      self$clock <- match.arg(clock)
      if (is.null(reset_states)){
        self$reset_states <- vector(mode = "double")
      } else{
        self$reset_states <- reset_states
      }
      
      # history
      self$start_state <- private$check_history(start_state)
      self$start_age <- private$check_history(start_age)
      
      # death state
      if (!is.null(death_state)){
        if (death_state > nrow(trans_mat)){
          stop("'death_state' cannot be larger than the number of rows in 'trans_mat'",
               call. = FALSE)
        } else{
          self$death_state <- death_state
        }
      } else{
        absorbing_states <- absorbing(trans_mat) 
        self$death_state <- absorbing_states[length(absorbing_states)]
      }
    },

    #' @description
    #' Simulate health state probabilities at distinct times by first 
    #' simulating trajectories through a multi-state model with 
    #' `IndivCtstm$sim_disease()`.
    #' @param t A numeric vector of times.
    #' @param ... Additional arguments to pass to `IndivCtstm$sim_disease()`.
    #' @return An object of class [stateprobs].  
    sim_stateprobs = function(t, ...){
      self$check()
      args <- c(trans_model = self, list(...))
      disprog <- do.call("indiv_ctstm_sim_disease", args)
      return(indiv_ctstm_sim_stateprobs(disprog, self, t))
    },
    
    #' @description
    #' Input validation for class. Checks that fields are the correct type. 
    check = function(){
      private$check_base()
    }
    
    
  ) # end public
  
)

# IndivCtstm -------------------------------------------------------------------
#' Individual-level continuous time state transition model
#'
#' @description
#' Simulate outcomes from an individual-level continuous time state transition 
#' model (CTSTM) from a fitted multi-state model. The class supports "clock-reset"
#' (i.e., semi-Markov), "clock-forward" (i.e., Markov), and mixtures of 
#' clock-reset and clock-forward models as described in 
#' [IndivCtstmTrans]. 
#' @format An [R6::R6Class] object.
#' @seealso [create_IndivCtstmTrans()], [IndivCtstmTrans]
#' @export
IndivCtstm <- R6::R6Class("IndivCtstm",
  private = list(  
    .stateprobs_ = NULL,
    .qalys_ = NULL,
    .costs_ = NULL,
    disprog_idx = NULL,
    
    sim_wlos = function(stateval_list, dr, stateval_type = c("costs", "qalys"),
                        sim_type, by_patient = FALSE, max_t = Inf,
                        lys = FALSE){
     
      stateval_type <- match.arg(stateval_type)
      if(is.null(self$disprog_)){
        stop("You must first simulate disease progression using '$sim_disease'.",
            call. = FALSE)
      }
      check_dr(dr)
      
      # Indexing patient and strategy ID's
      if (is.null(private$disprog_idx)){
        self$disprog_[, strategy_idx := .GRP, by = "strategy_id"]
        self$disprog_[, patient_idx := .GRP, by = "patient_id"]
        private$disprog_idx <- self$disprog_[, c("strategy_idx", "patient_idx"), with = FALSE]
        self$disprog_[, strategy_idx := NULL]
        self$disprog_[, patient_idx := NULL]
        private$disprog_idx[, strategy_idx := strategy_idx - 1]
        private$disprog_idx[, patient_idx := patient_idx - 1]
      }
      
      #  Categories
      n_cats <- length(stateval_list)
      if (stateval_type == "costs"){
        if (is.null(names(stateval_list))){
          categories <- paste0("Category ", seq(1, n_cats))
        } else{
          categories <- names(stateval_list)
        } # end if/else names for cost models
      } else{
        categories <- "qalys"
      } # end if/else costs vs. qalys
      
      
      # Maximum time
      if (!(length(max_t) %in% c(1, n_cats))){
        stop("'max_t' must either equal the number of cost categories or be of length 1.",
             call. = FALSE)
      }
      if (length(max_t) == 1){
        max_t <- rep(max_t, n_cats)
      }
      
      # Computation
      n_dr <- length(dr)
      wlos_list <- vector(mode = "list", length = n_cats * n_dr)
      counter <- 1
      for (i in 1:n_cats){
        for (j in 1:n_dr){
          C_wlos <- C_indiv_ctstm_wlos(self$disprog_, # Note: C++ re-indexing done at C level for disprog_
                                       private$disprog_idx$strategy_idx,
                                       private$disprog_idx$patient_idx,
                                       stateval_list[[i]], dr[j],
                                       sim_type, max_t[i])
          self$disprog_[, wlos := C_wlos]
          if (lys){
            C_los <- C_indiv_ctstm_los(self$disprog_, # Note: C++ re-indexing done at C level for disprog_
                                       private$disprog_idx$strategy_idx,
                                       private$disprog_idx$patient_idx,
                                       dr[j])
            self$disprog_[, lys := C_los]
            sdcols <- c("wlos", "lys")
          } else{
            sdcols <- "wlos"
          }
          if (by_patient == TRUE){
            by_cols <- c("sample", "strategy_id", "patient_id", "grp_id", "from")
            wlos_list[[counter]] <- self$disprog_[, lapply(.SD, sum), 
                                                        .SDcols = sdcols,
                                                        by = by_cols]
            setkeyv(wlos_list[[counter]], by_cols)
            # Pad missing health states within sample/strategy pairs with  NA's
            wlos_list[[counter]] <- wlos_list[[counter]][CJ(sample, strategy_id, patient_id, grp_id, from,
                                                              unique = TRUE)]
          } else{
            by_cols <- c("sample", "strategy_id", "grp_id", "from")
            n_patients <- self$trans_model$input_data$n_patients
            wlos_list[[counter]] <- self$disprog_[, lapply(.SD, sum), 
                                                        .SDcols = sdcols,
                                                        by = by_cols]
            wlos_list[[counter]][, wlos := wlos/n_patients]
            if(lys) wlos_list[[counter]][, lys := lys/n_patients]
            # Pad missing health states within sample/strategy pairs with NA's
            setkeyv(wlos_list[[counter]], by_cols)
            wlos_list[[counter]] <- wlos_list[[counter]][CJ(sample, strategy_id, grp_id, from,
                                                            unique = TRUE)]
          }
          wlos_list[[counter]][, "dr" := dr[j]]
          wlos_list[[counter]][, "category" := categories[i]]
          self$disprog_[, "wlos" := NULL]
          wlos_list[[counter]][, wlos := ifelse(is.na(wlos), 0, wlos)] # Replace padded NA's with 0's
          if (lys){
            self$disprog_[, "lys" := NULL]
            wlos_list[[counter]][, lys := ifelse(is.na(lys), 0, lys)] # Replace padded NA's with 0's
          }
          counter <- counter + 1
        }
      }
      wlos_dt <- rbindlist(wlos_list)
      setcolorder(wlos_dt, c(by_cols, "dr", "category", sdcols))
      setnames(wlos_dt, "from", "state_id")
      if (stateval_type == "costs"){
        setnames(wlos_dt, "wlos", "costs")
      } else{
        setnames(wlos_dt, "wlos", "qalys")
        wlos_dt[, category := NULL]
      }
      return(wlos_dt[,])
     } # end sim_wlos

  ), # end private  
                                                  
  public = list(
    #' @field trans_model The model for health state transitions. Must be an object 
    #' of class [IndivCtstmTrans]. 
    trans_model = NULL,
    
    #' @field utility_model The model for health state utility. Must be an object of
    #' class [StateVals].    
    utility_model = NULL,
    
    #' @field cost_models The models used to predict costs by health state. 
    #' Must be a list of objects of class [StateVals], where each element of the 
    #' list represents a different cost category.    
    cost_models = NULL,
    
    #' @field disprog_ A `data.table` simulated using `$sim_disease()` containing the following columns:
    #' * `sample`: A random sample from the PSA.
    #' * `strategy_id`: The treatment strategy ID.
    #' * `patient_id`: The patient ID.
    #' * `from`: The health state ID transitioned from.
    #' * `to`: The health state ID transitioned to.
    #' * `final`: An indicator equal to 1 if a patient is in their final health
    #' state during the simulation and 0 otherwise.
    #' * `time_start`: The time at the start of the interval.
    #' * `time_stop`: The time at the end of the interval.
    disprog_ = NULL,
    
    #' @field stateprobs_ An object of class [stateprobs] simulated using `$sim_stateprobs()`.    
    stateprobs_ = NULL,
    
    #' @field qalys_ An object of class [qalys] simulated using `$sim_qalys()`.
    qalys_ = NULL,
    
    #' @field costs_ An object of class [costs] simulated using `$sim_costs()`.
    costs_ = NULL,
    
    #' @description
    #' Create a new `IndivCtstm` object.
    #' @param trans_model The `trans_model` field.
    #' @param utility_model The `utility_model` field.
    #' @param cost_models The `cost_models` field.
    #' @return A new `IndivCtstm` object.      
    initialize = function(trans_model = NULL, utility_model = NULL, cost_models = NULL) {
      self$trans_model <- trans_model
      self$utility_model = utility_model
      self$cost_models = cost_models
    },
    
    #' @description
    #' Simulate disease progression (i.e., individual trajectories through a multi-state
    #' model using an individual patient simulation).
    #' @param max_t  A scalar or vector denoting the length of time to simulate the model. 
    #' If a vector, must be equal to the number of simulated patients. 
    #' @param max_age A scalar or vector denoting the maximum age to simulate each patient until.
    #'  If a vector, must be equal to the number of simulated patients.
    #' @param progress An integer, specifying the PSA iteration (i.e., sample) that should be 
    #' printed every `progress` PSA iterations. For example, if `progress = 2`, 
    #' then every second PSA iteration is printed. Default is `NULL`, 
    #' in which case no output is printed.
    #' @return An instance of self with simulated output stored in `disprog_`.
    sim_disease = function(max_t = 100, max_age = 100, progress = NULL){
      if(!inherits(self$trans_model, "IndivCtstmTrans")){
        stop("'trans_model' must be an object of class 'IndivCtstmTrans'",
          call. = FALSE)
      }
      self$trans_model$check()
      
      self$qalys_ <- NULL
      self$costs_ <- NULL
      self$stateprobs_ <- NULL
      private$disprog_idx <- NULL
      self$disprog_ <- indiv_ctstm_sim_disease(self$trans_model,
                                               max_t = max_t,
                                               max_age = max_age,
                                               progress = progress)
      self$stateprobs_ <- NULL
      invisible(self)
    },
    
    #' @description
    #' Simulate health state probabilities as a function of time using the 
    #' simulation output stored in `disprog`.
    #' @param t TA numeric vector of times.
    #' @return An instance of `self` with simulated output of class [stateprobs] 
    #' stored in `stateprobs_`.    
    sim_stateprobs = function(t){
      if(is.null(self$disprog_)){
        stop("You must first simulate disease progression using '$sim_disease'.",
            call. = FALSE)
      }
      self$stateprobs_ <- indiv_ctstm_sim_stateprobs(self$disprog_,
                                                     self$trans_model,
                                                     t = t)
      invisible(self)
    },
    
    #' @description
    #' Simulate quality-adjusted life-years (QALYs) using [sim_qalys()].
    #' @param dr Discount rate.
    #' @param type `"predict"` for mean values or `"random"` for random samples 
    #' as in `$sim()` in [StateVals].
    #' @param lys If `TRUE`, then life-years are simulated in addition to QALYs.
    #' @param by_patient If `TRUE`, then QALYs are computed at the patient level.
    #'  If `FALSE`, then QALYs are averaged across patients by health state.
    #' @return An instance of `self` with simulated output stored
    #' in `qalys_`.    
    sim_qalys = function(dr = .03, type = c("predict", "random"),
                         lys = TRUE,
                         by_patient = FALSE){
      if(!inherits(self$utility_model, "StateVals")){
        stop("'utility_model' must be an object of class 'StateVals'",
          call. = FALSE)
      }
      type <- match.arg(type)
      self$qalys_ <- private$sim_wlos(list(self$utility_model), dr, "qalys", type, by_patient,
                                          lys = lys)
      invisible(self)
    },
    
    #' @description
    #' Simulate costs using [sim_costs()].
    #' @param dr Discount rate.
    #' @param type `"predict"` for mean values or `"random"` for random samples 
    #' as in `$sim()` in [StateVals].
    #' @param by_patient If `TRUE`, then QALYs are computed at the patient level.
    #'  If `FALSE`, then QALYs are averaged across patients by health state.
    #' @param max_t  Maximum time duration to compute costs once a patient has
    #'   entered a (new) health state. By default, equal to `Inf`, 
    #'   so that costs are computed over the entire duration that a patient is in
    #'   a given health state. If time varies by each cost category, then time 
    #'   can also be passed as a numeric vector of length equal to the number of 
    #'   cost categories (e.g., `c(1, 2, Inf, 3)` for a model with four cost categories).
    #' @return An instance of `self` with simulated output stored
    #' in `costs_`.    
    sim_costs = function(dr = .03, type = c("predict", "random"), by_patient = FALSE, max_t = Inf){
      if(!is.list(self$cost_models)){
        stop("'cost_models' must be a list of objects of class 'StateVals'",
          call. = FALSE)
      } else{
        for (i in 1:length(self$cost_models)){
          if (!inherits(self$cost_models[[i]], "StateVals")){
            stop("'cost_models' must be a list of objects of class 'StateVals'",
                call. = FALSE)
          }
        }
      }
      type <- match.arg(type)
      self$costs_ <- private$sim_wlos(self$cost_models, dr, "costs", type, by_patient, max_t)
      invisible(self)
    },
    
    #' @description
    #' Summarize costs and QALYs so that cost-effectiveness analysis can be performed. 
    #' See [summarize_ce()].    
    summarize = function() {
      check_summarize(self)
      summarize_ce(self$costs_, self$qalys_)
    }
    
  ) # end public
) # end class
