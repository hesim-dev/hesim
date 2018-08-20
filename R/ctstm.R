# CtstmTrans -------------------------------------------------------------------
CtstmTrans <- R6::R6Class("CtstmTrans",
  private = list(
    check_base = function(){
      if(!inherits(self$data, "input_data")){
        stop("'data' must be an object of class 'input_data'",
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
        if(nrow(self$trans_mat) != length(self$params)){
          stop(paste0("The number of models in 'params' must equal the number of rows/columns ",
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
    data = NULL,
    params = NULL,
    trans_mat = NULL,
    
    hazard = function(t){
      private$summary(t, "hazard")
    },
    
    cumhazard = function(t){
       private$summary(t, "cumhazard")
    }
  )
)

# IndivCtstmTrans --------------------------------------------------------------
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
  disprog$sim[, strategy_index := .GRP, by = "strategy_id"]
  strategy_index <- disprog$sim$strategy_id - 1
  disprog$sim[, strategy_index := NULL]

  ## Dimensions of simulation
  if (!"line" %in% names(disprog$sim)){
    n_lines <- 1
  } # to do after incorporating treatment lines: case where there are multiple treatment lines
  n_strategies <- length(disprog$unique_strategy_id)
  n_patients <- length(disprog$unique_patient_id)
  
  ## Computation
  stprobs <- C_ctstm_indiv_stateprobs(disprog$sim, t, 
                                      disprog$n_samples,
                                      n_strategies, 
                                      disprog$unique_strategy_id,
                                      strategy_index,
                                      disprog$n_states,
                                      n_patients,
                                      n_lines)
  stprobs <- data.table(stprobs)
      
  ## C++ to R indexing
  stprobs[, "sample" := get("sample") + 1]
  stprobs[, "state_id" := get("state_id") + 1]
  
  # Return
  return(stprobs[])      
}

#' Create \code{IndivCtstmTrans} object
#' 
#' \code{mTrans} is a generic function for creating an object of class
#' \code{\link{IndivCtstmTrans}}.
#' @param object A fitted statistical model. 
#' @param data An object of class "expanded_hesim_data" returned by 
#' \code{\link{expand_hesim_data}}.
#' @param n Number of random observations of the parameters to draw.
#' @param trans_mat The transition matrix describing the states and transitions in a 
#' multi-state model in the format from the \link[mstate]{mstate} package. See \link{IndivCtstmTrans}.
#' @param point_estimate If \code{TRUE}, then the point estimates are returned and and no samples are drawn.
#' @param ... Further arguments passed to \code{IndivCtstmTrans$new()} in \code{\link{IndivCtstmTrans}}.
#' @return Returns an \code{\link{R6Class}} object of class \code{\link{IndivCtstmTrans}}.
#' @seealso \code{\link{IndivCtstmTrans}}
#' @name create_IndivCtstmTrans
#' @rdname create_IndivCtstmTrans
#' @export
create_IndivCtstmTrans <- function(object, data, trans_mat, n = 1000, point_estimate = FALSE, ...){
  UseMethod("create_IndivCtstmTrans", object)
}

#' @export
#' @rdname create_IndivCtstmTrans
create_IndivCtstmTrans.flexsurvreg_list <- function(object, data, trans_mat, n = 1000, point_estimate = FALSE, ...){
  input_data <- create_input_data(object, data, id_vars = c("strategy_id", "patient_id"))
  params <- create_params(object, n = n, point_estimate = point_estimate)
  return(IndivCtstmTrans$new(data = input_data, params = params, trans_mat = trans_mat, ...))
}

#' @export
#' @rdname create_IndivCtstmTrans
create_IndivCtstmTrans.flexsurvreg <- function(object, data, trans_mat, n = 1000, point_estimate = FALSE, ...){
  input_data <- create_input_data(object, data, id_vars = c("strategy_id", "patient_id", "transition_id"))
  params <- create_params(object, n = n, point_estimate = point_estimate)
  return(IndivCtstmTrans$new(data = input_data, params = params, trans_mat = trans_mat, ...))
}

#' @export
IndivCtstmTrans <- R6::R6Class("IndivCtstmTrans",
  inherit = CtstmTrans,
  
  private = list(
    .death_state = NULL
  ), # end private
  
  active = list(
    death_state = function(value) {
      if (missing(value)) {
        private$.death_state
      } else {
        stop("'$death_state' is read only", call. = FALSE)
      }
     }
   ), # end active  
  
  public = list(
    start_state = NULL,
    start_time = NULL,
    start_age = NULL,
    
    
    initialize = function(data, params, trans_mat, 
                          start_state = 1,
                          start_time = 0,
                          start_age = 38,
                          death_state = NULL) {
      self$data <- data
      self$params <- params
      self$trans_mat <- trans_mat
      
      # history
      if (length(start_age) == 1){
        self$start_age <- rep(start_age, data$n_patients)
      }
      if (length(self$start_age) != data$n_patients){
        stop("The length of 'start_age' must equal 'n_patients' in 'data'.",
             call. = FALSE)
      } 
      self$start_state <- start_state
      
      
      # death state
      if (!is.null(death_state)){
        if (death_state > nrow(trans_mat)){
          stop("'death_state' cannot be larger than the number of rows in 'trans_mat'",
               call. = FALSE)
        } else{
          private$.death_state <- death_state
        }
      } else{
        absorbing_states <- absorbing(trans_mat) 
        private$.death_state <- absorbing_states[length(absorbing_states)]
      }
    },
    
    sim_disease = function(max_t = 100, max_age = 100){
      private$check_base()
      
      # Simulate
      disprog <- C_ctstm_sim_disease(self, self$start_state, self$start_age,
                                     self$death_state - 1, max_t, max_age)
      disprog <- data.table(disprog)
      disprog[, sample := sample + 1]
      disprog[, from := from + 1]
      disprog[, to := to + 1]
      disprog[, line := NULL] # to do after incorporating treatment lines: case where there are multiple treatment lines
      
      if (inherits(self$params, "params_surv_list")){
        n_samples <- self$params[[1]]$n_samples
      } else{
        n_samples <- self$params$n_samples
      }      
      out <- list(sim = disprog[,], 
                  n_samples = n_samples,
                  n_states = nrow(self$trans_mat),
                  unique_strategy_id = unique(self$data$strategy_id),
                  unique_patient_id = unique(self$data$patient_id))
      class(out) <- "indiv_ctstm_disprog"
      return(out)
    },
    
    
    sim_stateprobs = function(disprog = NULL, t, ...){
      return(indiv_ctstm_sim_stateprobs(disprog, self, t))
    },
    
    check = function(){
      private$check_base()
    }
    
    
  ) # end public
  
)

# IndivCtstm -------------------------------------------------------------------
#' @export
IndivCtstm <- R6::R6Class("IndivCtstm",
  private = list(  
    .disprog_ = NULL,
    .stateprobs_ = NULL,
    .qalys_ = NULL,
    .costs_ = NULL,
    disprog_idx = NULL,
    
    sim_wlos = function(stateval_list, dr, stateval_type = c("costs", "qalys"),
                        sim_type){
     
      stateval_type <- match.arg(stateval_type)
      if(is.null(self$disprog_)){
        stop("You must first simulate disease progression using '$sim_disease'.",
            call. = FALSE)
      }      
      
      # Indexing patient and strategy ID's
      if (is.null(private$disprog_idx)){
        self$disprog_$sim[, strategy_idx := .GRP, by = "strategy_id"]
        self$disprog_$sim[, patient_idx := .GRP, by = "patient_id"]
        private$disprog_idx <- self$disprog_$sim[, c("strategy_idx", "patient_idx"), with = FALSE]
        self$disprog_$sim[, strategy_idx := NULL]
        self$disprog_$sim[, patient_idx := NULL]
        private$disprog_idx[, strategy_idx := strategy_idx - 1]
        private$disprog_idx[, patient_idx := patient_idx - 1]
      }
      
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
      

      n_dr <- length(dr)
      wlos_list <- vector(mode = "list", length = n_cats * n_dr)
      counter <- 1
      for (i in 1:n_cats){
        for (j in 1:n_dr){
          C_wlos <- C_indiv_ctstm_wlos(self$disprog_$sim, # Note: C++ re-indexing done at C level for disprog_
                                       private$disprog_idx$strategy_idx,
                                       private$disprog_idx$patient_idx,
                                       stateval_list[[i]], dr[j],
                                       sim_type)
          self$disprog_$sim[, wlos := C_wlos]
          wlos_list[[counter]] <- self$disprog_$sim[, .(wlos = sum(wlos)), 
                                       by = c("sample", "strategy_id", "patient_id")]
          wlos_list[[counter]][, "dr" := dr[j]]
          wlos_list[[counter]][, "category" := categories[i]]
          self$disprog_$sim[, "wlos" := NULL]
          counter <- counter + 1
        }
      }
      wlos_dt <- rbindlist(wlos_list)
      setcolorder(wlos_dt, c("sample", "strategy_id", "patient_id", "dr", "category", "wlos"))
      if (stateval_type == "costs"){
        setnames(wlos_dt, "wlos", "costs")
      } else{
        setnames(wlos_dt, "wlos", "qalys")
        wlos_dt[, category := NULL]
      }
      return(wlos_dt[,])
     } # end sim_wlos

  ), # end private  
                                                  
  active = list(
    disprog_ = function(value) {
      if (missing(value)) {
        private$.disprog_
      } else {
        stop("'$disprog_' is read only", call. = FALSE)
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
  ), # end active
  
  public = list(
    trans_model = NULL,
    utility_model = NULL,
    cost_models = NULL,
    initialize = function(trans_model = NULL, disprog = NULL, utility_model = NULL, cost_models = NULL) {
      self$trans_model <- trans_model
      self$utility_model = utility_model
      self$cost_models = cost_models
      private$.disprog_ = disprog
    },
    
    sim_disease = function(max_t = 100, max_age = 100){
      if(!inherits(self$trans_model, "IndivCtstmTrans")){
        stop("'trans_model' must be an object of class 'IndivCtstmTrans'",
          call. = FALSE)
      }
      self$trans_model$check()
      
      private$.qalys_ <- NULL
      private$.costs_ <- NULL
      private$.stateprobs_ <- NULL
      private$disprog_idx <- NULL
      private$.disprog_ <- self$trans_model$sim_disease(max_t = max_t,
                                                        max_age = max_age)
      private$.stateprobs_ <- NULL
      invisible(self)
    },
    
    sim_stateprobs = function(t){
      if(is.null(self$disprog_)){
        stop("You must first simulate disease progression using '$sim_disease'.",
            call. = FALSE)
      }
      private$.stateprobs_ <- indiv_ctstm_sim_stateprobs(self$disprog_, t = t)
      invisible(self)
    },
    
    sim_qalys = function(dr = .03, type = c("predict", "random")){
      if(!inherits(self$utility_model, "StateVals")){
        stop("'utility_model' must be an object of class 'StateVals'",
          call. = FALSE)
      }
      type <- match.arg(type)
      private$.qalys_ <- private$sim_wlos(list(self$utility_model), dr, "qalys", type)
      invisible(self)
    },
    
    sim_costs = function(dr = .03, type = c("predict", "random")){
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
      private$.costs_ <- private$sim_wlos(self$cost_models, dr, "costs", type)
      invisible(self)
    },
    
    summarize = function(stat = mean) {
      if (is.null(self$costs_)) {
        stop("Cannot summarize costs without first simulating 'costs_' with '$sim_costs()'.",
               call. = FALSE)
      }
      
      if (is.null(self$qalys_)) {
        stop("Cannot summarize QALYs without first simulating 'qalys_' with '$sim_qalys()'.",
              call. = FALSE)
      }      
      
      costs_summary <- self$costs_[, lapply(.SD, stat), by = c("category", "dr", "sample", "strategy_id"),
                                   .SDcols = "costs"]
      qalys_summary <- self$qalys_[, lapply(.SD, stat), by = c("dr", "sample", "strategy_id"),
                                   .SDcols = "qalys"]
      ce <- list(costs = costs_summary, qalys = qalys_summary)
      class(ce) <- "ce"
      return(ce)
    }
    
  ) # end public
) # end class
