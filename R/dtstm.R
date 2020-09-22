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
    .start_stateprobs = NULL,
    .trans_mat = NULL,
    
    set_start_stateprobs = function(x){
      if (any(x < 0)){
        stop("All elements of 'state_stateprobs' must be non-negative.")
      }
      if (any(is.infinite(x))){
        stop("Elements of 'state_stateprobs' cannot be infinite.")
      }
      if (all(x == 0)) {
        x_len <- length(x)
        private$.start_stateprobs <- rep(1/x_len, x_len)
      } else{
        private$.start_stateprobs <- x/sum(x)
      }
    },
    
    set_trans_mat = function(x){
      if (is.null(x)){
        private$.trans_mat <- x
      } else{
        if (!is.matrix(x)){
          stop("'trans_mat' must be a matrix.")
        }
        K <- ncol(x)
        for (i in 1:K) {
          NA_count <- sum(is.na(x[i, ]))
          max_value <- ifelse(NA_count, NA, K - 1 - NA_count)
          if (!is.na(max_value)){
            if (!all(x[i, ] == 0:max_value)){
              stop(paste0("'trans_mat' is not of the correct form. Each row should ",
                          "contain integers from 0 to K - 1 where K is the number of ",
                          "possible transitions (i.e., non-NA elements)"))
            }
          } 
        }
        private$.trans_mat <- x
      }
    },
    
    get_n_states = function(){
      if (is.null(self$input_data)){
        return(ncol(self$params$value[,,1]))
      } else{
        return(nrow(self$trans_mat))
      }
    }
  ),        
  
  active = list(
    #' @field start_stateprobs A non-negative vector with length equal to the number of
    #' health states containing the probability that the cohort is in each health
    #' state at the start of the simulation. For example, 
    #' if there were three states and the cohort began the simulation in state 1,
    #' then `start_stateprobs = c(1, 0, 0)`. Automatically normalized to sum to 1.
    #' If `NULL`, then a vector with the first element equal to 1 and 
    #' all remaining elements equal to 0.
    start_stateprobs = function(x) {
      if (missing(x)) return(private$.start_stateprobs)
      private$set_start_stateprobs(x)
    },
    
    #' @field trans_mat A transition matrix describing the states and transitions
    #' in a discrete-time multi-state model. Only required if the model is 
    #' parameterized using multinomial logistic regression. The `(i,j)` element 
    #' represents a transition from state `i` to state `j`. Each possible transition 
    #' from row `i` should be based on a separate multinomial logistic regression
    #' and ordered from `0` to `K - 1` where `K` is the number of 
    #' possible transitions. Transitions that are not possible should be `NA`.
    #' and the reference category for each row should be `0`.    
    trans_mat = function(x){
      if(missing(x)) return(private$.trans_mat)
      private$set_trans_mat(x)
    }
  ),

  public = list(
    #' @field params Parameters for simulating health state transitions.
    #' Supports objects of class [tparams_transprobs] or [params_mlogit].
    params = NULL,
    
    #' @field input_data An object of class [input_mats].
    input_data = NULL,
    
    #' @field cycle_length The length of a model cycle in terms of years.
    #' The default is `1` meaning that model cycles are 1 year long.
    cycle_length = NULL,
    
    #' @description
    #' Create a new `CohortDtstmTrans` object.
    #' @param params The `params` field.
    #' @param input_data The `input_data` field.
    #' @param trans_mat The `trans_mat` field.
    #' @param start_stateprobs The `start_stateprobs` field.
    #' @param cycle_length The `cycle_length` field.
    #' @return A new `CohortDtstmTrans` object.
    initialize = function(params,
                          input_data = NULL,
                          trans_mat = NULL,
                          start_stateprobs = NULL, 
                          cycle_length = 1){
      self$params <- params
      self$input_data <- input_data
      private$set_trans_mat(trans_mat)
      if (!is.null(input_data) & is.null(trans_mat)){
        stop("If 'input_data' is not NULL, then 'trans_mat' cannot be NULL")
      }
      if (!is.null(start_stateprobs)){
        private$set_start_stateprobs(start_stateprobs)
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
      stprobs <- C_cohort_dtstm_sim_stateprobs(self, times)
      stprobs <- data.table(stprobs)
      stprobs[, sample := sample + 1]
      stprobs[, state_id := state_id + 1]
      check_patient_wt(self, stprobs)
      setattr(stprobs, "class", 
              c("stateprobs", "data.table", "data.frame"))
      return(stprobs[])
    }
  )
)

#' Create \code{CohortDtstmTrans} object
#' 
#' A generic function for creating an object of class `CohortDtstmTrans`.
#' @param object An object of the appropriate class. 
#' @param ... Further arguments passed to `CohortDtstmTrans$new()` in 
#' \code{CohortDtstmTrans}.
#' @param input_data An object of class `expanded_hesim_data` returned by 
#' [expand.hesim_data()]
#' @param trans_mat A transition matrix describing the states and transitions 
#' in a discrete-time multi-state model. See [CohortDtstmTrans].
#' @param n Number of random observations of the parameters to draw.
#' @param point_estimate If `TRUE`, then the point estimates are returned and and no samples are drawn.
#' @export
create_CohortDtstmTrans <- function(object, ...){
  UseMethod("create_CohortDtstmTrans", object)
} 

#' @export
#' @rdname create_CohortDtstmTrans
create_CohortDtstmTrans.multinom_list <- function(object, input_data,
                                                  trans_mat,
                                                  n = 1000, point_estimate = FALSE,
                                                  ...){
  input_mats <- create_input_mats(object, input_data)
  params <- create_params(object, n = n, point_estimate = point_estimate)
  return(CohortDtstmTrans$new(params = params, input_data = input_mats, 
                              trans_mat = trans_mat, ...))
}

# CohortDtstm ------------------------------------------------------------------
#' Cohort discrete time state transition model
#'
#' @description
#' Simulate outcomes from a cohort discrete time state transition model.
#' @examples 
#' library("data.table")
#' library("nnet")
#' transitions_data <- data.table(multinom3_exdata$transitions)
#'
#' # Treatment strategies, target population, and model structure
#' n_patients <- 4
#'patients <- transitions_data[year == 1, .(patient_id, age, female)][
#'  sort(sample.int(nrow(transitions_data[year == 1]), n_patients))][
#'    , grp_id := 1:n_patients]
#' hesim_dat <- hesim_data(
#'   patients = patients,
#'   strategies = data.table(strategy_id = 1:2,
#'                           strategy_name = c("Reference", "Intervention")),
#'   states = data.table(state_id = c(1, 2),
#'                       state_name = c("Healthy", "Sick")) # Non-death health states
#' )
#' tmat <- rbind(c(0, 1, 2),
#'               c(NA, 0, 1),
#'               c(NA, NA, NA))
#'
#' # Parameter estimation
#' ## Multinomial logistic regression
#' data_healthy <- transitions_data[state_from == "Healthy"]
#' fit_healthy <- multinom(state_to ~ strategy_name + female + age,
#'                         data = data_healthy, trace = FALSE)
#' data_sick <- droplevels(transitions_data[state_from == "Sick"])
#' fit_sick <- multinom(state_to ~ strategy_name + female + age,
#'                      data = data_sick, trace = FALSE)
#' transfits <- multinom_list(healthy = fit_healthy, sick = fit_sick)
#'
#' ## Utility
#' utility_tbl <- stateval_tbl(multinom3_exdata$utility,
#'                             dist = "beta",
#'                             hesim_data = hesim_dat)
#'
#' ## Costs
#' drugcost_tbl <- stateval_tbl(multinom3_exdata$costs$drugs,
#'                              dist = "fixed",
#'                              hesim_data = hesim_dat) 
#' medcost_tbl <- stateval_tbl(multinom3_exdata$costs$medical,
#'                             dist = "gamma",
#'                             hesim_data = hesim_dat)  
#'
#'# Economic model
#' n_samples <- 3
#'
#' ## Construct model
#' ### Transitions
#' transmod_data <- expand(hesim_dat)
#' transmod <- create_CohortDtstmTrans(transfits,
#'                                     input_data = transmod_data,
#'                                     trans_mat = tmat,
#'                                     n = n_samples,
#'                                     point_estimate = FALSE)
#'
#' ### Utility
#' utilitymod <- create_StateVals(utility_tbl, n = n_samples)
#'
#' ### Costs
#' drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples)
#' medcostmod <- create_StateVals(medcost_tbl, n = n_samples)
#' costmods <- list(Drug = drugcostmod,
#'                  Medical = medcostmod)
#'
#' ### Combine
#' econmod <- CohortDtstm$new(trans_model = transmod,
#'                            utility_model = utilitymod,
#'                            cost_models = costmods)
#'
#' ## Simulate outcomes
#' econmod$sim_stateprobs(n_cycles = 20)
#' econmod$sim_qalys(dr = .03)
#' econmod$sim_costs(dr = .03)
#' econmod$summarize()
#' econmod$summarize(by_grp = TRUE)
#'
#' @format An [R6::R6Class] object.
#' @seealso [create_CohortDtstm()], [CohortDtstmTrans], 
#' [create_CohortDtstmTrans()]
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
    #' Simulate health state probabilities using `CohortDtstmTrans$sim_stateprobs()`.
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
    #' `utility_model`. See `vignette("expected-values")` for details.
    #' @param dr Discount rate.
    #' @param integrate_method Method used to integrate state values when computing (QALYs).
    #' @param lys If `TRUE`, then life-years are simulated in addition to QALYs.
    #' @return An instance of `self` with simulated output of class [qalys] stored
    #' in `qalys_`.
    sim_qalys = function(dr = .03,
                         integrate_method = c("trapz", "riemann_left", "riemann_right"),
                         lys = TRUE){
      self$qalys_ <- sim_qalys(self$stateprobs_, self$utility_model, dr, integrate_method,
                               lys)
      invisible(self)
    },
    
    #' @description
    #' Simulate costs as a function of `stateprobs_` and `cost_models`. 
    #' See `vignette("expected-values")` for details.
    #' @param dr Discount rate.
    #' @param integrate_method Method used to integrate state values when computing costs.
    #' @return An instance of `self` with simulated output of class [costs] stored
    #' in `costs_`.
    sim_costs = function(dr = .03, 
                         integrate_method = c("trapz", "riemann_left", "riemann_right")){
      self$costs_ <- sim_costs(self$stateprobs_, self$cost_models, dr, integrate_method)
      invisible(self)
    },

    #' @description
    #' Summarize costs and QALYs so that cost-effectiveness analysis can be performed. 
    #' See [summarize_ce()].    
    #' @param by_grp If `TRUE`, then costs and QALYs are computed by subgroup. If
    #' `FALSE`, then costs and QALYs are aggregated across all patients (and subgroups).    
    summarize = function(by_grp = FALSE) {
      check_summarize(self)
      return(summarize_ce(self$costs_, self$qalys_, by_grp))
    }
  )
)

#' Create \code{CohortDtstm} object
#' 
#' A generic function for creating an object of class [CohortDtstm].
#' @param object An object of the appropriate class. 
#' @param input_data 	An object of class [expanded_hesim_data][expand.hesim_data()].
#' @param cost_args A list of further arguments passed to `StateVals$new()` in 
#' [StateVals] when initiating cost models.
#' @param utility_args A list of further arguments passed to `StateVals$new()` in
#' [StateVals] when initiating the utility model.
#' @param ... Further arguments passed to `CohortDtstmTrans$new()` in 
#' [CohortDtstmTrans]. 
#'  
#' @export
create_CohortDtstm <- function(object, ...){
  UseMethod("create_CohortDtstm", object)
} 

#' @export
#' @rdname create_CohortDtstm
create_CohortDtstm.model_def <- function(object, input_data, cost_args = NULL,
                                         utility_args = NULL, ...){
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
    utility_model <- create_StateVals(model_inputs, 
                                      cost = FALSE,
                                      init_args = utility_args)
  }
  
  ## Cost models
  n_cost_models <- length(model_inputs$cost)
  if (n_cost_models > 0){
   cost_models <- vector(mode = "list", length = length(model_inputs$cost))
   names(cost_models) <- names(model_inputs$cost)
   for (i in 1:n_cost_models){
     cost_models[[i]] <- create_StateVals(model_inputs, 
                                          name = names(cost_models)[i],
                                          init_args = cost_args[[i]]) 
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

