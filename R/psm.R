# PsmCurves --------------------------------------------------------------------

#' Create \code{PsmCurves} object
#' 
#' \code{create_PsmCurves} is a function for creating an object of class
#' \code{\link{PsmCurves}}.
#' @param object Fitted survival models.
#' @param input_data An object of class `expanded_hesim_data` returned by 
#' \code{\link{expand.hesim_data}}. Must be expanded by the data tables "strategies" and
#' "patients". 
#' @param n Number of random observations of the parameters to draw.
#' @param uncertainty Method determining how parameter uncertainty should be handled. See
#'  documentation in [`create_params()`].
#' @param est_data A `data.table` or `data.frame` of estimation data 
#' used to fit survival models during bootstrap replications.
#' @param ... Further arguments passed to or from other methods. Passed to [`create_params.partsurvfit()`]
#' when `object` is of class [`flexsurvreg_list`].
#' @return Returns an `R6Class` object of class [`PsmCurves`].
#' @seealso [`PsmCurves`]
#' @export
create_PsmCurves <- function(object, ...){
  UseMethod("create_PsmCurves", object)
} 
 
#' @export
#' @rdname create_PsmCurves
create_PsmCurves.flexsurvreg_list <- function(object, input_data, n = 1000, 
                                              uncertainty = c("normal", "bootstrap", "none"), 
                                              est_data = NULL,
                                               ...){
  # For backwards compatibility until deprecated point_estimate argument is no longer supported
  is_uncertainty_missing <- missing(uncertainty)
  uncertainty <- deprecate_point_estimate(list(...)$point_estimate, uncertainty,
                                          is_uncertainty_missing)
  uncertainty <- deprecate_bootstrap(list(...)$bootstrap, uncertainty,
                                     is_uncertainty_missing)
  
  # Code to always keep
  uncertainty <- match.arg(uncertainty)
  if (uncertainty == "bootstrap" & is.null(est_data)){
    stop("If uncertainty == 'bootstrap', then est_data cannot be NULL")
  }
  psfit <- partsurvfit(object, est_data)
  input_mats <- create_input_mats(psfit, input_data, id_vars = c("strategy_id", "patient_id"))
  params <- create_params(psfit, n = n, uncertainty = uncertainty)
  return(PsmCurves$new(input_data = input_mats, params = params))
}

#' @export
#' @rdname create_PsmCurves
create_PsmCurves.params_surv_list <- function(object, input_data, ...){
  input_mats <- create_input_mats(object, input_data)
  return(PsmCurves$new(input_data = input_mats, params = object))
}


#' Partitioned survival curves
#'
#' @description
#' Summarize N-1 survival curves for an N-state partitioned survival model.
#' @format An [R6::R6Class] object.
#' @examples 
#' library("flexsurv")
#'
#' # Simulation data
#' dt_strategies <- data.frame(strategy_id = c(1, 2, 3))
#' dt_patients <- data.frame(patient_id = seq(1, 3),
#'                           age = c(45, 50, 60),
#'                           female = c(0, 0, 1))
#' hesim_dat <- hesim_data(strategies = dt_strategies,
#'                         patients = dt_patients)
#'
#' # Fit survival models
#' surv_est_data <- psm4_exdata$survival
#' fit1 <- flexsurv::flexsurvreg(Surv(endpoint1_time, endpoint1_status) ~ age,
#'                               data = surv_est_data, dist = "exp")
#' fit2 <- flexsurv::flexsurvreg(Surv(endpoint2_time, endpoint2_status) ~ age,
#'                               data = surv_est_data, dist = "exp")
#' fit3 <- flexsurv::flexsurvreg(Surv(endpoint3_time, endpoint3_status) ~ age,
#'                               data = surv_est_data, dist = "exp")
#' fits <- flexsurvreg_list(fit1, fit2, fit3)
#'
#' # Form PsmCurves
#' surv_input_data <- expand(hesim_dat, by = c("strategies", "patients"))
#' psm_curves <- create_PsmCurves(fits, input_data = surv_input_data, n = 3,
#'                                uncertainty = "bootstrap", est_data = surv_est_data)
#'
#' # Summarize survival curves
#' head(psm_curves$quantile(p = c(.25, .5, .75)))
#' head(psm_curves$survival(t = seq(0, 3, by = .1)))
#' head(psm_curves$rmst(t = c(2, 5)))
#' 
#' @seealso [`Psm`], [`create_PsmCurves()`]
#' @export
PsmCurves <- R6::R6Class("PsmCurves",
  private = list(
    summary = function(x, type = c("hazard", "cumhazard", "survival", 
                                   "rmst", "quantile"), 
                       dr = 0){
      self$check()
      type <- match.arg(type)
      res <- data.table(C_psm_curves_summary(self, x, type, dr))
      res[, curve := curve + 1]
      res[, sample := sample + 1]
      check_patient_wt(self, res)
      if (type %in% c("hazard", "cumhazard", "survival", "rmst")){
        setnames(res, "x", "t")
      } else if (type == "quantile"){
        setnames(res, "x", "p")
      }
      if (type == "hazard") setnames(res, "value", "hazard")
      if (type == "cumhazard") setnames(res, "value", "cumhazard")
      if (type == "survival") setnames(res, "value", "survival")
      if (type == "rmst") setnames(res, "value", "rmst")
      if (type == "quantile") setnames(res, "value", "quantile")
      return(res[])
    }
  ),                            
                              
  public = list(
    #' @field params An object of class [`params_surv_list`].
    params = NULL,
    
    #' @field input_data An object of class [`input_mats`]. Each row in `X` must
    #' be a unique treatment strategy and patient.
    input_data = NULL,

    #' @description
    #' Create a new `PsmCurves` object.
    #' @param params The `params` field.
    #' @param input_data The `input_data` field.
    #' @return A new `PsmCurves` object.
    initialize = function(params, input_data) {
      self$params <- params
      self$input_data <- input_data
    },

    #' @description
    #' Predict the hazard function for each survival curve as a function of time.
    #' @param t  A numeric vector of times.
    #' @return A `data.table` with columns `sample`, `strategy_id`,
    #' `patient_id`, `grp_id`, `curve` (the curve number), `t`, and `hazard`.
    hazard = function(t){
      return(private$summary(x = t, type = "hazard"))
    },
    
    #' @description
    #' Predict the cumulative hazard function for each survival curve as a function of time.
    #' @param t  A numeric vector of times.
    #' @return A `data.table` with columns `sample`, `strategy_id`,
    #' `patient_id`, `grp_id`, `curve`, `t`, and `cumhazard`.
    cumhazard = function(t){
      return(private$summary(x = t, type = "cumhazard"))
    },
    
    #' @description
    #' Predict survival probabilities for each survival curve as a function of time.
    #' @param t  A numeric vector of times.
    #' @return An object of class [`survival`].   
    survival = function(t){
      res <- private$summary(x = t, type = "survival")
      setattr(res, "class", 
              c("survival", "data.table", "data.frame"))
      return(res)
    },

    #' @description
    #' Predict the restricted mean survival time up until time points `t`
    #'  for each survival curve.
    #' @param t  A numeric vector of times.
    #' @param dr Discount rate.
    #' @return A `data.table` with columns `sample`, `strategy_id`,
    #' `patient_id`, `grp_id`, `curve`, `t`, and `rmst`.     
    rmst = function(t, dr = 0){
      return(private$summary(x = t, type = "rmst", dr = dr))
    },
    
    #' @description
    #' Predict quantiles of the survival distribution for each survival curve.
    #' @param p  A numeric vector of probabilities for computing quantiles.
    #' @return A `data.table` with columns `sample`, `strategy_id`,
    #' `patient_id`, `grp_id`, `curve`, `p` and `quantile`.       
    quantile = function(p){
      return(private$summary(x = p, type = "quantile"))
    },

    #' @description
    #' Input validation for class. Checks that fields are the correct type.     
    check = function(){
      if(!inherits(self$input_data, "input_mats")){
        stop("'input_data' must be an object of class 'input_mats'",
            call. = FALSE)
      }
      if(!inherits(self$params, c("params_surv_list"))){
        stop("Class of 'params' is not supported. See documentation.",
             call. = FALSE)
      }
    }
  )
)

# Psm --------------------------------------------------------------------------
#' N-state partitioned survival model
#'
#' @description
#' Simulate outcomes from an N-state partitioned survival model.
#' @format An [R6::R6Class] object.
#' @param dr Discount rate.
#' @param integrate_method Method used to integrate state values when computing 
#' costs or QALYs. Options are `trapz` for the trapezoid rule,
#' `riemann_left` for a left Riemann sum, and  
#' `riemann_right` for a right Riemann sum.
#' @examples
#' library("flexsurv")
#' library("ggplot2")
#' theme_set(theme_bw())
#'
#' # Simulation data
#' strategies <- data.frame(strategy_id = c(1, 2, 3), 
#'                          strategy_name = paste0("Strategy ", 1:3))
#' patients <- data.frame(patient_id = seq(1, 3),
#'                        age = c(45, 50, 60),
#'                        female = c(0, 0, 1))
#' states <- data.frame(state_id =  seq(1, 3),
#'                      state_name = paste0("State ", seq(1, 3)))
#' hesim_dat <- hesim_data(strategies = strategies,
#'                         patients = patients,
#'                         states = states)
#' labs <- c(
#'   get_labels(hesim_dat),
#'   list(curve = c("Endpoint 1" = 1,
#'                 "Endpoint 2" = 2,
#'                  "Endpoint 3" = 3))
#' )
#' n_samples <- 2
#'
#' # Survival models
#' surv_est_data <- psm4_exdata$survival
#' fit1 <- flexsurv::flexsurvreg(Surv(endpoint1_time, endpoint1_status) ~ factor(strategy_id),
#'                               data = surv_est_data, dist = "exp")
#' fit2 <- flexsurv::flexsurvreg(Surv(endpoint2_time, endpoint2_status) ~ factor(strategy_id),
#'                               data = surv_est_data, dist = "exp")
#' fit3 <- flexsurv::flexsurvreg(Surv(endpoint3_time, endpoint3_status) ~ factor(strategy_id),
#'                               data = surv_est_data, dist = "exp")
#' fits <- flexsurvreg_list(fit1, fit2, fit3)
#'
#' surv_input_data <- expand(hesim_dat, by = c("strategies", "patients"))
#' psm_curves <- create_PsmCurves(fits, input_data = surv_input_data,
#'                                uncertainty = "bootstrap", est_data = surv_est_data, 
#'                                n = n_samples)
#'
#' # Cost model(s)
#' cost_input_data <- expand(hesim_dat, by = c("strategies", "patients", "states"))
#' fit_costs_medical <- stats::lm(costs ~ female + state_name, 
#'                                data = psm4_exdata$costs$medical)
#' psm_costs_medical <- create_StateVals(fit_costs_medical, 
#'                                       input_data = cost_input_data, 
#'                                       n = n_samples)
#'
#' # Utility model
#' utility_tbl <- stateval_tbl(tbl = data.frame(state_id = states$state_id,
#'                                              min = psm4_exdata$utility$lower,
#'                                              max = psm4_exdata$utility$upper),
#'                             dist = "unif")
#' psm_utility <- create_StateVals(utility_tbl, n = n_samples,
#'                                 hesim_data = hesim_dat)
#'
#' # Partitioned survival decision model
#' psm <- Psm$new(survival_models = psm_curves,
#'                utility_model = psm_utility,
#'                cost_models = list(medical = psm_costs_medical))
#' psm$sim_survival(t = seq(0, 5, 1/12))
#' autoplot(psm$survival_, labels = labs, ci = FALSE, ci_style = "ribbon")
#' psm$sim_stateprobs()
#' autoplot(psm$stateprobs_, labels = labs)
#' psm$sim_costs(dr = .03)
#' head(psm$costs_)
#' head(psm$sim_qalys(dr = .03)$qalys_)
#'
#' @seealso [`PsmCurves`], [`create_PsmCurves()`]
#' @references [Incerti and Jansen (2021)](https://arxiv.org/abs/2102.09437).
#' See Section 2.3 for a mathematical description of a PSM and Section 4.2 for an 
#' example in oncology. The mathematical approach used to simulate costs and QALYs from
#' state probabilities is described in Section 2.1.
#' @export
Psm <- R6::R6Class("Psm",
  public = list(
    #' @field survival_models The survival models used to predict survival curves. Must be
    #' an object of class [`PsmCurves`].
    survival_models = NULL,
    
    #' @field utility_model The model for health state utility. Must be an object of
    #' class [`StateVals`].
    utility_model = NULL,
    
    #' @field cost_models The models used to predict costs by health state. 
    #' Must be a list of objects of class [`StateVals`], where each element of the 
    #' list represents a different cost category.    
    cost_models = NULL,
    
    #' @field n_states Number of states in the partitioned survival model.
    n_states = NULL,
    
    #' @field t_ A numeric vector of times at which survival curves were predicted. Determined
    #' by the argument `t` in `$sim_curves()`.
    t_ = NULL,
    
    #' @field survival_ An object of class [survival] simulated using `sim_survival()`.
    survival_ = NULL,
    
    #' @field stateprobs_ An object of class [stateprobs] simulated using `$sim_stateprobs()`.
    stateprobs_ = NULL,
    
    #' @field qalys_ An object of class [qalys] simulated using `$sim_qalys()`.
    qalys_ = NULL,
    
    #' @field costs_ An object of class [costs] simulated using `$sim_costs()`.
    costs_ = NULL,

    #' @description
    #' Create a new `Psm` object.
    #' @param survival_models The `survival_models` field.
    #' @param utility_model The `utility_model` field.
    #' @param cost_models The `cost_models` field.
    #' @details `n_states` is set equal to the number of survival models plus one.
    #' @return A new `Psm` object.    
    initialize = function(survival_models = NULL, utility_model = NULL, cost_models = NULL) {
      self$survival_models <- survival_models
      self$cost_models = cost_models
      self$utility_model = utility_model
      self$n_states <- length(self$survival_models$params) + 1
    },
    
    #' @description
    #' Simulate survival curves as a function of time using `PsmCurves$survival()`.
    #' @param t A numeric vector of times. The first element must be `0`.
    #' @return An instance of `self` with simulated output from `PsmCurves$survival()`
    #' stored in `survival_`.
    sim_survival = function(t){
      if (t[1] !=0){
        stop("The first element of 't' must be 0.", call. = FALSE)
      }
      if(!inherits(self$survival_models, "PsmCurves")){
        stop("'survival_models' must be of class 'PsmCurves'.")
      }
      self$survival_models$check()
      self$survival_ <- self$survival_models$survival(t)
      setattr(self$survival_, "class", 
              c("survival", "data.table", "data.frame"))
      setattr(self$survival_, "size", 
              c(get_size(self$survival_models),
                n_states = self$n_states,
                n_times = length(t)))
      self$t_ <- t
      self$stateprobs_ <- NULL
      invisible(self)
    },
    
    #' @description
    #' Simulate health state probabilities from `survival_` using a partitioned
    #' survival analysis.
    #' @return An instance of `self` with simulated output of class [stateprobs] 
    #' stored in `stateprobs_`.    
    sim_stateprobs = function(){
      if(is.null(self$survival_)){
        stop("You must first simulate survival curves using '$sim_survival'.",
            call. = FALSE)
      }
      self$stateprobs_ <- sim_stateprobs(self$survival_)
      invisible(self)
    },

    #' @description
    #' Simulate quality-adjusted life-years (QALYs) as a function of `stateprobs_` and
    #' `utility_model`. 
    #' @param lys If `TRUE`, then life-years are simulated in addition to QALYs.
    #' @return An instance of `self` with simulated output of class [qalys] stored
    #' in `qalys_`.    
    sim_qalys = function(dr = .03, integrate_method = c("trapz", "riemann_left", "riemann_right"),
                         lys = TRUE){
      self$qalys_ <- sim_qalys(self$stateprobs_, self$utility_model, dr, 
                               integrate_method, lys)
      invisible(self)
    },
    
    #' @description
    #' Simulate costs as a function of `stateprobs_` and `cost_models`. 
    #' See `vignette("expected-values")` for details.
    #' @return An instance of `self` with simulated output of class [costs] stored
    #' in `costs_`.    
    sim_costs = function(dr = .03, integrate_method = c("trapz", "riemann_left", "riemann_right")){
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