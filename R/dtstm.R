# CohortDtstmTrans -------------------------------------------------------------
#' Transitions for a cohort discrete time state transition model
#'
#' @description
#' Simulate health state transitions in a cohort discrete time state transition model.
#' @format An [R6::R6Class] object.
#' @seealso [create_CohortDtstmTrans()], [CohortDtstm]
#' @export
CohortDtstmTrans <- R6::R6Class("CohortDtstmTrans",
  private = list(
    get_n_states = function(){
      if (is.null(self$input_mats)){
        return(ncol(self$params$value[,,1]))
      }
    }
  ),                             

  public = list(
    #' @field params Parameters for simulating health state transitions. Currently
    #' supports objects of class [tparams_transprobs].
    params = NULL,
    
    #' @field input_mats An object of class [input_mats].
    input_mats = NULL,
    
    #' @field start_stateprobs A vector with length equal to the number of
    #' health states containing the probability that the cohort is in each health
    #' state at the start of the simulation. For example, 
    #' if there were three states and the cohort began the simulation in state 1,
    #' then `start_stateprobs = c(1, 0, 0)`. If `NULL``, then a vector with the
    #' first element equal to 1 and all remaining elements equal to 0.
    start_stateprobs = NULL,
    
    #' @field cycle_length The length of a model cycle in terms of years.
    #' The default is `1` meaning that model cycles are 1 year long.
    cycle_length = NULL,
    
    #' @description
    #' Create a new `CohortDtstmTrans` object.
    #' @param params The `params` field.
    #' @param input_mats The `input_mats` field.
    #' @param start_stateprobs The `start_stateprobs` field.
    #' @param cycle_length The `cycle_length` field.
    #' @return A new `CohortDtstmTrans` object.
    initialize = function(params,
                          input_mats = NULL,
                          start_stateprobs = NULL, 
                          cycle_length = 1){
      self$params <- params
      self$input_mats <- input_mats
      if (!is.null(start_stateprobs)){
        self$start_stateprobs <- start_stateprobs
      } else{
        self$start_stateprobs <- c(1, rep(0, private$get_n_states() - 1))
      }
      self$cycle_length <- cycle_length
    },
    
    #' @description
    #' Simulate probability of being in each health state during each model cycle. 
    #' @param n_cycles The number of model cycles to simulate the model for.
    #' @return An object of class [stateprobs].
    sim_stateprobs = function(n_cycles){
      times <- seq(0, n_cycles/self$cycle_length, length.out = n_cycles + 1)
      stprobs <- C_cohort_dtstm_sim_stateprobs(self, 
                                               times,
                                               self$params$n_samples)
      stprobs <- data.table(stprobs)
      stprobs[, sample := sample + 1]
      stprobs[, state_id := state_id + 1]
      check_patient_wt(self, stprobs)
      return(stprobs[])
    }
  )
)

#' Create \code{CohortDtstmTrans} object
#' 
#' A generic function for creating an object of class \code{CohortDtstmTrans}.
#' @param object An object of the appropriate class. 
#' @param ... Further arguments passed to \code{CohortDtstmTrans$new()} in 
#' \code{CohortDtstmTrans}.
#'  
#' @export
create_CohortDtstmTrans <- function(object, ...){
  UseMethod("create_CohortDtstmTrans", object)
} 

# CohortDtstm ------------------------------------------------------------------
#' Cohort discrete time state transition model
#'
#' @description
#' Simulate outcomes from a cohort discrete time state transition model.
#' @format An [R6::R6Class] object.
#' @seealso [CohortDtstmTrans], [create_CohortDtstmTrans()]
#' @export
CohortDtstm <- R6::R6Class("CohortDtstm",
  public = list(
    #' @field trans_model The model for health state transitions. Must be an object 
    #' of class [CohortDtstmTrans]. 
    trans_model = NULL,
    
    #' @field utility_model The model for health state utility. Must be an object of
    #' class [StateVals].
    utility_model = NULL,
    
    #' @field cost_models The models used to predict costs by health state. 
    #' Must be a list of objects of class [StateVals], where each element of the 
    #' list represents a different cost category.
    cost_models = NULL,
    
    #' @field stateprobs_ An object of class [stateprobs] simulated using `$sim_stateprobs()`.
    stateprobs_ = NULL,
    
    #' @field qalys_ An object of class [qalys] simulated using `$sim_qalys()`.
    qalys_ = NULL,
    
    #' @field costs_ An object of class [costs] simulated using `$sim_costs()`.
    costs_ = NULL,
    
    #' @description
    #' Create a new `CohortDtstm` object.
    #' @param trans_model The `trans_model` field.
    #' @param utility_model The `utility_model` field.
    #' @param cost_models The `cost_models` field.
    #' @return A new `CohortDtstm` object.
    initialize = function(trans_model = NULL, utility_model = NULL, cost_models = NULL) {
      self$trans_model <- trans_model
      self$utility_model = utility_model
      self$cost_models = cost_models
    },
    
    #' @description
    #' Simulate health state probabilities using `CohortCtstmTrans$sim_stateprobs()`.
    #' @param n_cycles The number of model cycles to simulate the model for.
    #' @return An instance of `self` with simulated output of class [stateprobs] 
    #' stored in `stateprobs_`.
    sim_stateprobs = function(n_cycles){
      self$stateprobs_ <- self$trans_model$sim_stateprobs(n_cycles)
      setattr(self$stateprobs_, "class", 
              c("stateprobs", "data.table", "data.frame"))
      invisible(self)
    },
    
    #' @description
    #' Simulate quality-adjusted life-years (QALYs) as a function of `stateprobs_` and
    #' `utility_model`. See 
    #' the [vignette](https://hesim-dev.github.io/hesim/dev/articles/wlos.html) for details.
    #' @param dr Discount rate.
    #' @param method Method used to integrate state values when computing (QALYs).
    #' @param lys If `TRUE`, then life-years are simulated in addition to QALYs.
    #' @return An instance of `self` with simulated output of class [qalys] stored
    #' in `qalys_`.
    sim_qalys = function(dr = .03,
                         method = c("trapz", "riemann_left", "riemann_right"),
                         lys = FALSE){
      self$qalys_ <- sim_qalys(self$stateprobs_, self$utility_model, dr, method,
                               lys)
      invisible(self)
    },
    
    #' @description
    #' Simulate costs as a function of `stateprobs_` and `cost_models`. 
    #' See the [vignette](https://hesim-dev.github.io/hesim/dev/articles/wlos.html) for details.
    #' @param dr Discount rate.
    #' @param method Method used to integrate state values when computing costs.
    #' @return An instance of `self` with simulated output of class [costs] stored
    #' in `costs_`.
    sim_costs = function(dr = .03, 
                         method = c("trapz", "riemann_left", "riemann_right")){
      self$costs_ <- sim_costs(self$stateprobs_, self$cost_models, dr, method)
      invisible(self)
    },

    #' @description
    #' Summarize costs and QALYs so that cost-effectiveness analysis can be performed. 
    #' See [summarize_ce()].        
    summarize = function() {
      check_summarize(self)
      return(summarize_ce(self$costs_, self$qalys_))
    }
  )
)

#' Create \code{CohortDtstm} object
#' 
#' A generic function for creating an object of class [CohortDtstm].
#' @param object An object of the appropriate class. 
#' @param input_data 	An object of class [expanded_hesim_data][expand.hesim_data()].
#' @param ... Further arguments passed to `CohortDtstmTrans$new()` in 
#' [CohortDtstmTrans]. 
#'  
#' @export
create_CohortDtstm <- function(object, ...){
  UseMethod("create_CohortDtstm", object)
} 

#' @export
#' @rdname create_CohortDtstm
create_CohortDtstm.model_def <- function(object, input_data, ...){
  model_inputs <- eval_model(object, input_data)
  
  # Individual models
  ## Transition model
  if (is.null(model_inputs$tpmatrix)){
    trans_model <- NULL
  } else{
    tparams <- tparams_transprobs(model_inputs)
    trans_model <- CohortDtstmTrans$new(params = tparams, ...) 
  }
  
  ## Utility model
  if (is.null(model_inputs$utility)){
    utility_model <- NULL
  } else{
    utility_model <- create_StateVals(model_inputs, cost = FALSE) 
  }
  
  ## Cost models
  n_cost_models <- length(model_inputs$cost)
  if (n_cost_models > 0){
   cost_models <- vector(mode = "list", length = length(model_inputs$cost))
   names(cost_models) <- names(model_inputs$cost)
   for (i in 1:n_cost_models){
     cost_models[[i]] <- create_StateVals(model_inputs, name = names(cost_models)[i])
   }
  } else{
    cost_models <- NULL
  }
  
  # Combine
  econ_model <- CohortDtstm$new(trans_model = trans_model,
                                utility_model = utility_model,
                                cost_models = cost_models)
  
  
  return(econ_model)
}
