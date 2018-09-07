# PsmCurves --------------------------------------------------------------------

#' Create \code{PsmCurves} object
#' 
#' \code{create_PsmCurves} is a function for creating an object of class
#' \code{\link{PsmCurves}} from an object of class \code{\link{partsurvfit}}.
#' @param object An object of class \code{\link{partsurvfit}}.
#' @param data An object of class "expanded_hesim_data" returned by 
#' \code{\link{expand.hesim_data}}. Must be expanded by the data tables "strategies" and
#' "patients". 
#' @param n Number of random observations of the parameters to draw.
#' @param point_estimate If \code{TRUE}, then the point estimates are returned and and no samples are drawn.
#' @param bootstrap If TRUE, then \code{n} bootstrap replications are drawn by refitting the survival
#'  models in \code{object} on resamples of the sample data; if FALSE, then the parameters for each survival
#'  model are independently draw from multivariate normal distributions.  
#' @param ... Further arguments passed to or from other methods. Currently unused. 
#' @return Returns an \code{\link{R6Class}} object of class \code{\link{PsmCurves}}.
#' @seealso \code{\link{PsmCurves}}
#' @export
create_PsmCurves <- function(object, data, n = 1000, point_estimate = FALSE,
                                bootstrap = TRUE){
  if (!inherits(object, c("partsurvfit"))){
    stop("'Object' must be of class 'partsurvfit'.")
  }
  input_data <- create_input_data(object, data, id_vars = c("strategy_id", "patient_id"))
  params <- create_params(object, n = n, point_estimate = point_estimate, bootstrap = bootstrap)
  return(PsmCurves$new(data = input_data, params = params))
}

# Manual documentation in PsmCurves.Rd
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
    data = NULL,
    params = NULL,

    initialize = function(data, params) {
      self$data <- data
      self$params <- params
    },
    
    hazard = function(t){
      return(private$summary(x = t, type = "hazard"))
    },
    
    cumhazard = function(t){
      return(private$summary(x = t, type = "cumhazard"))
    },
    
    survival = function(t){
      return(private$summary(x = t, type = "survival"))
    },
    
    rmst = function(t, dr = 0){
      return(private$summary(x = t, type = "rmst", dr = dr))
    },
    
    quantile = function(p){
      return(private$summary(x = p, type = "quantile"))
    },
    
    check = function(){
      if(!inherits(self$data, "input_data")){
        stop("'data' must be an object of class 'input_data'",
            call. = FALSE)
      }
      if(!inherits(self$params, c("params_surv_list", 
                                  "joined_params_surv_list"))){
        stop("Class of 'params' is not supported. See documentation.",
             call. = FALSE)
      }
    }
  )
)

# Psm --------------------------------------------------------------------------
# Manual documentation in Psm.Rd
#' @export
Psm <- R6::R6Class("Psm",
  private = list(
    .t_ = NULL,
    .survival_ = NULL,
    .stateprobs_ = NULL,
    .costs_ = NULL,
    .qalys_ = NULL,
    
    sim_wlos = function(dr, type){
      if(is.null(self$stateprobs_)){
        stop("You must first simulate health state probabilities using '$sim_stateprobs'.",
             call. = FALSE)
      }
      
      statvalmods <- switch(type,
                           costs = self$cost_models,
                           qalys = list(self$utility_model))
      
      statvalmods_name <- switch(type,
                                costs = "cost_models",
                                qalys = "utility_model")
      
      # Check number of samples
      expected_samples <- max(self$stateprobs_$sample)
      for (i in 1:length(statvalmods)){
        if (statvalmods[[i]]$params$n_samples != expected_samples){
          msg <- paste0("Number of samples in '", statvalmods_name, "' must equal to ",
                        " the number of samples in 'survival_models', which is ",
                         expected_samples)
          stop(msg, call. = FALSE)
        }
      }
      
      # Check number of states
      for (i in 1:length(statvalmods)){
        if(self$n_states != statvalmods[[i]]$data$n_states + 1){
          msg <- paste0("The number of survival models must equal the number of states in '",
                        statvalmods_name, "' - 1.")
          stop(msg, call. = FALSE)
        }
      } # loop over models
      
      stateprobs <- self$stateprobs_[state_id != self$n_states] 
      
      if (type == "costs"){
        if (is.null(names(self$cost_models))){
          categories <- paste0("Type ", seq(1, length(self$cost_models)))
        } else{
            categories <- names(self$cost_models)
        } # end if/else names for cost models
      } else{
        categories <- "qalys"
      } # end if/else costs vs. qalys

      res <- data.table(C_psm_sim_wlos(self, stateprobs, dr, type, categories))
      res[, sample := sample + 1]
      return(res[])
    } # end sim_wlos()
  ), # end private
  
  active = list(
    t_ = function(value) {
      if (missing(value)) {
        private$.t_
      } else {
        stop("'$t_' is read only", call. = FALSE)
      }
    },
    
    survival_ = function(value) {
      if (missing(value)) {
        private$.survival_
      } else {
        stop("'$survival_' is read only", call. = FALSE)
      }
    },
    
    stateprobs_ = function(value) {
      if (missing(value)) {
        private$.stateprobs_
      } else {
        stop("'$stateprobs_' is read only", call. = FALSE)
      }
    },
    
    qalys_ = function(value) {
      if (missing(value)) {
        private$.qalys_
      } else {
        stop("'$qalys_' is read only", call. = FALSE)
      }
    },
    
    costs_ = function(value) {
      if (missing(value)) {
        private$.costs_
      } else {
        stop("'$costs_' is read only", call. = FALSE)
      }
    }
    
  ),
                        
  public = list(
    survival_models = NULL,
    utility_model = NULL,
    cost_models = NULL,
    n_states = NULL,

    initialize = function(survival_models, utility_model = NULL, cost_models = NULL) {
      self$survival_models <- survival_models
      self$cost_models = cost_models
      self$utility_model = utility_model
      self$n_states <- length(self$survival_models$params) + 1
    },
    
    sim_survival = function(t){
      if (t[1] !=0){
        stop("The first element of 't' must be 0.", call. = FALSE)
      }
      if(!inherits(self$survival_models, "PsmCurves")){
        stop("'survival_models' must be of class 'PsmCurves'.")
      }
      self$survival_models$check()
      private$.survival_ <- self$survival_models$survival(t)
      private$.t_ <- t
      private$.stateprobs_ <- NULL
      invisible(self)
    },
    
    sim_stateprobs = function(){
      if(is.null(self$survival_)){
        stop("You must first simulate survival curves using '$sim_survival'.",
            call. = FALSE)
      }
      res <- C_psm_sim_stateprobs(self$survival_,
                                  n_samples = self$survival_models$params[[1]]$n_samples,
                                  n_strategies = self$survival_models$data$n_strategies,
                                  n_patients = self$survival_models$data$n_patients,
                                  n_states = self$n_states,
                                  n_times = length(self$t_))
      prop_cross <- res$n_crossings/nrow(res$stateprobs)
      if (prop_cross > 0){
        warning(paste0("Survival curves crossed ", round(prop_cross * 100, 2), 
                       " percent of the time."),
                call. = FALSE)
      }
      stateprobs <- data.table(res$stateprobs)
      stateprobs[, state_id := state_id + 1]
      stateprobs[, sample := sample + 1]
      private$.stateprobs_ <- stateprobs[]
      invisible(self)
    },
    
    sim_qalys = function(dr = .03){
      self$utility_model$check()
      qalys <- private$sim_wlos(dr, type = "qalys")
      setnames(qalys, "value", "qalys")
      private$.qalys_ <- qalys
      invisible(self)
    },
    
    sim_costs = function(dr = .03){
      if(!is.list(self$cost_models)){
        stop("'cost_models' must be a list", call. = FALSE)
      }
      for (i in 1:length(self$cost_models)){
        self$cost_models[[i]]$check()
      }
      costs <- private$sim_wlos(dr, type = "costs")
      setnames(costs, "value", "costs")
      private$.costs_ <- costs
      invisible(self)
    }
  )
)