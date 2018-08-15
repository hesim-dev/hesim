# CtstmTrans -------------------------------------------------------------------
#' Create \code{CtstmTrans} object
#' 
#' \code{create_CtstmTrans} is a generic function for creating an object of class
#' \code{\link{CtstmTrans}}.
#' @param object A fitted statistical model. 
#' @param data An object of class "expanded_hesim_data" returned by 
#' \code{\link{expand_hesim_data}}.
#' @param n Number of random observations of the parameters to draw.
#' @param trans_mat The transition matrix describing the states and transitions in a 
#' multi-state model in the format from the \link[mstate]{mstate} package. See \link{CtstmTrans}.
#' @param point_estimate If \code{TRUE}, then the point estimates are returned and and no samples are drawn.
#' @param ... Further arguments passed to \code{CtstmTrans$new()} in \code{\link{CtstmTrans}}.
#' @return Returns an \code{\link{R6Class}} object of class \code{\link{CtstmTrans}}.
#' @seealso \code{\link{CtstmTrans}}
#' @name create_CtstmTrans
#' @rdname create_CtstmTrans
#' @export
create_CtstmTrans <- function(object, data, trans_mat, n = 1000, point_estimate = FALSE, ...){
  UseMethod("create_CtstmTrans", object)
}

#' @export
#' @rdname create_CtstmTrans
create_CtstmTrans.flexsurvreg_list <- function(object, data, trans_mat, n = 1000, point_estimate = FALSE, ...){
  input_data <- create_input_data(object, data, id_vars = c("strategy_id", "patient_id"))
  params <- create_params(object, n = n, point_estimate = point_estimate)
  return(CtstmTrans$new(data = input_data, params = params, trans_mat = trans_mat, ...))
}

#' @export
#' @rdname create_CtstmTrans
create_CtstmTrans.flexsurvreg <- function(object, data, trans_mat, n = 1000, point_estimate = FALSE, ...){
  input_data <- create_input_data(object, data, id_vars = c("strategy_id", "patient_id", "transition_id"))
  params <- create_params(object, n = n, point_estimate = point_estimate)
  return(CtstmTrans$new(data = input_data, params = params, trans_mat = trans_mat, ...))
}

indiv_ctstm_stateprobs <- function(disease_prog, t, trans_model) {
  if (inherits(trans_model$params, "params_surv_list")){
    n_samples <- trans_model$params[[1]]$n_samples
  } else{
    n_samples <- trans_model$params$n_samples
  }
  if (is.null(trans_model$data$n_lines)){
    n_lines <- 1
  } # to do after incorporating treatment lines: case where there are multiple treatment lines
  
  stprobs <- C_ctstm_indiv_stateprobs(disease_prog, t, n_samples,
                                      trans_model$data$n_strategies, 
                                      nrow(trans_model$trans_mat),
                                      trans_model$data$n_patients,
                                      n_lines)
  
  stprobs <- data.table(stprobs)
      
  ## C++ to R indexing
  stprobs[, "sample" := get("sample") + 1]
  stprobs[, "strategy_id" := get("strategy_id") + 1]
  stprobs[, "state_id" := get("state_id") + 1]
      
  # Return
  return(stprobs[])
}

#' @export
CtstmTrans <- R6::R6Class("CtstmTrans",
  private = list(
      .death_state = NULL,
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
    data = NULL,
    params = NULL,
    trans_mat = NULL,
    start_ages = NULL,
    initialize = function(data, params, trans_mat, start_ages = rep(38, data$n_patients),
                          death_state = NULL) {
      self$data <- data
      self$params <- params
      self$trans_mat <- trans_mat
      
      # starting ages
      if (length(start_ages) != data$n_patients){
        stop("The length of 'start_ages' must equal 'n_patients' in 'data'.",
             call. = FALSE)
      } else{
       self$start_ages <- start_ages 
      }
      
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
    
    hazard = function(t){
      private$summary(t, "hazard")
    },
    
    cumhazard = function(t){
       private$summary(t, "cumhazard")
    },
    
    sim_stateprobs = function(t, start_state = 1, max_t = 100, max_age = 100,
                              clock = "reset"){
      self$check()
      disprog <- C_ctstm_sim_disease(self, start_state, self$start_ages,
                                     self$death_state - 1, max_t, max_age)
      return(indiv_ctstm_stateprobs(disprog, t, self))
    },
    
    check = function(){
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
    }
  )
)

# IndivCtstm -------------------------------------------------------------------
#' @export
IndivCtstm <- R6::R6Class("IndivCtstm",
  private = list(  
    .disease_prog_ = NULL,
    .stateprobs_ = NULL,
    .qalys_ = NULL,
    .costs_ = NULL,
    
    create_C_disease_prog = function(R_disease_prog){
      x <- copy(R_disease_prog)
      x[, sample := sample - 1]
      x[, strategy_id := strategy_id - 1]
      x[, patient_id := patient_id - 1]
      x[, from := from - 1]
      x[, to := to - 1]
      return(x)
    },
  
    sim_wlos = function(stateval_list, dr, stateval_type = c("costs", "qalys"),
                        sim_type){
      C_disease_prog <- private$create_C_disease_prog(self$disease_prog_)
      stateval_type <- match.arg(stateval_type)
      if(is.null(self$disease_prog_)){
        stop("You must first simulate disease progression using '$sim_disease'.",
            call. = FALSE)
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
          C_wlos <- C_indiv_ctstm_wlos(C_disease_prog, stateval_list[[i]], dr[j],
                                       sim_type)
          self$disease_prog_[, wlos := C_wlos]
          wlos_list[[counter]] <- self$disease_prog_[, .(wlos = sum(wlos)), 
                                       by = c("sample", "strategy_id", "patient_id")]
          wlos_list[[counter]][, "dr" := dr[j]]
          wlos_list[[counter]][, "category" := categories[i]]
          self$disease_prog_[, "wlos" := NULL]
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
    disease_prog_ = function(value) {
      if (missing(value)) {
        private$.disease_prog_
      } else {
        stop("'$disease_prog_' is read only", call. = FALSE)
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
    initialize = function(trans_model, utility_model = NULL, cost_models = NULL) {
      self$trans_model <- trans_model
      self$utility_model = utility_model
      self$cost_models = cost_models
    },
    
    sim_disease = function(max_t = 100, max_age = 100){
      private$.qalys_ <- NULL
      private$.costs_ <- NULL
      private$.stateprobs_ <- NULL
      start_state <- 1
      self$trans_model$check()
      sim <- C_ctstm_sim_disease(self$trans_model, start_state, 
                                 self$trans_model$start_ages, 
                                 self$trans_model$death_state - 1, max_t, max_age)
      sim <- data.table(sim)
      
      # C++ to R indexing
      sim[, sample := sample + 1]
      sim[, strategy_id := strategy_id + 1]
      sim[, line := NULL]
      sim[, patient_id := patient_id + 1]
      sim[, from := from + 1]
      sim[, to := to + 1]
      
      # Update class states
      private$.disease_prog_ <- sim[]
      private$.stateprobs_ <- NULL
      invisible(self)
    },
    
    sim_stateprobs = function(t){
      if(is.null(self$disease_prog_)){
        stop("You must first simulate disease progression using '$sim_disease'.",
            call. = FALSE)
      }
      
      C_disease_prog <- private$create_C_disease_prog(self$disease_prog_) # convert to C++ indices
      private$.stateprobs_ <- indiv_ctstm_stateprobs(C_disease_prog, t, self$trans_model)
      invisible(self)
    },
    
    sim_qalys = function(dr = .03, type = c("predict", "random")){
      type <- match.arg(type)
      private$.qalys_ <- private$sim_wlos(list(self$utility_model), dr, "qalys", type)
      invisible(self)
    },
    
    sim_costs = function(dr = .03, type = c("predict", "random")){
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
