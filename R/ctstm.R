# CtstmTrans -------------------------------------------------------------------
CtstmTrans <- R6::R6Class("CtstmTrans",
  private = list(
    check_base = function(){
      if(!inherits(self$data, "input_mats")){
        stop("'data' must be an object of class 'input_mats'",
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
indiv_ctstm_sim_disease <- function(trans_model, max_t = 100, max_age = 100){
  sample <- from <- to <- line <- NULL # to avoid no visible bindings CRAN warning
  
  # Simulate
  disprog <- C_ctstm_sim_disease(trans_model, trans_model$start_state - 1, 
                                 trans_model$start_age,
                                 trans_model$start_time,
                                 trans_model$death_state - 1, 
                                 max_t, max_age)
  disprog <- data.table(disprog)
  disprog[, sample := sample + 1]
  disprog[, from := from + 1]
  disprog[, to := to + 1]
  disprog[, line := NULL] # to do after incorporating treatment lines: case where there are multiple treatment lines
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
  strategy_index <- disprog$strategy_id - 1
  disprog[, strategy_index := NULL]

  ## Dimensions of simulation
  if (!"line" %in% names(disprog)){
    n_lines <- 1
  } # to do after incorporating treatment lines: case where there are multiple treatment lines
  n_states <- nrow(trans_model$trans_mat)
  n_strategies <- trans_model$data$n_strategies
  n_patients <- trans_model$data$n_patients
  if (inherits(trans_model$params, "params_surv_list")){
    n_samples <- trans_model$params[[1]]$n_samples
  } else{
    n_samples <- trans_model$params$n_samples
  }
  unique_strategy_id <- unique(trans_model$data$strategy_id)
  
  ## Computation
  stprobs <- C_ctstm_indiv_stateprobs(disprog, 
                                      t, 
                                      n_samples,
                                      n_strategies, 
                                      unique_strategy_id,
                                      strategy_index,
                                      n_states,
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
#' A generic function for creating an object of class \code{\link{IndivCtstmTrans}}.
#' @param object A fitted statistical model. 
#' @param data An object of class "expanded_hesim_data" returned by 
#' \code{\link{expand.hesim_data}}.
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
  X <- create_input_mats(object, data)
  params <- create_params(object, n = n, point_estimate = point_estimate)
  return(IndivCtstmTrans$new(data = X, params = params, trans_mat = trans_mat, ...))
}

#' @export
#' @rdname create_IndivCtstmTrans
create_IndivCtstmTrans.flexsurvreg <- function(object, data, trans_mat, n = 1000, point_estimate = FALSE, ...){
  X <- create_input_mats(object, data)
  params <- create_params(object, n = n, point_estimate = point_estimate)
  return(IndivCtstmTrans$new(data = X, params = params, trans_mat = trans_mat, ...))
}

#' @export
#' @rdname create_IndivCtstmTrans
create_IndivCtstmTrans.params_surv <- function(object, data, trans_mat, ...){
  X <- create_input_mats(object, data)
  return(IndivCtstmTrans$new(data = X, params = object, trans_mat = trans_mat, ...))
}

#' @export
IndivCtstmTrans <- R6::R6Class("IndivCtstmTrans",
  inherit = CtstmTrans,
  
  private = list(
    .death_state = NULL,
    
    check_history = function(field){
      field_name <- deparse(substitute(field))
      if (length(field) !=1 & length(field) != self$data$n_patients){
        stop(paste0("The length of '", field_name, "' must either be 1 or the number ",
                          "of simulated patients."),
                 call. = FALSE)        
      }
      if (length(field) == 1){
        field <- rep(field, self$data$n_patients)
      }
      return(field)
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
      self$start_state <- private$check_history(start_state)
      self$start_time <- private$check_history(start_time)
      self$start_age <- private$check_history(start_age)
      
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
    
    sim_stateprobs = function(t, ...){
      self$check()
      args <- c(trans_model = self, list(...))
      disprog <- do.call("indiv_ctstm_sim_disease", args)
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
                        sim_type, by_patient = FALSE, max_t = Inf,
                        lys = FALSE){
     
      stateval_type <- match.arg(stateval_type)
      if(is.null(self$disprog_)){
        stop("You must first simulate disease progression using '$sim_disease'.",
            call. = FALSE)
      }      
      
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
            by_cols <- c("sample", "strategy_id", "patient_id", "from")
            wlos_list[[counter]] <- self$disprog_[, lapply(.SD, sum), 
                                                        .SDcols = sdcols,
                                                        by = by_cols]
            setkeyv(wlos_list[[counter]], by_cols)
            # Pad missing health states within sample/strategy pairs with  NA's
            wlos_list[[counter]] <- wlos_list[[counter]][CJ(sample, strategy_id, patient_id, from,
                                                              unique = TRUE)]
          } else{
            by_cols <- c("sample", "strategy_id", "from")
            n_patients <- self$trans_model$data$n_patients
            wlos_list[[counter]] <- self$disprog_[, lapply(.SD, sum), 
                                                        .SDcols = sdcols,
                                                        by = by_cols]
            wlos_list[[counter]][, wlos := wlos/n_patients]
            if(lys) wlos_list[[counter]][, lys := lys/n_patients]
            # Pad missing health states within sample/strategy pairs with NA's
            setkeyv(wlos_list[[counter]], by_cols)
            wlos_list[[counter]] <- wlos_list[[counter]][CJ(sample, strategy_id, from,
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
    initialize = function(trans_model = NULL, utility_model = NULL, cost_models = NULL) {
      self$trans_model <- trans_model
      self$utility_model = utility_model
      self$cost_models = cost_models
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
      private$.disprog_ <- indiv_ctstm_sim_disease(self$trans_model,
                                                   max_t = max_t,
                                                   max_age = max_age)
      private$.stateprobs_ <- NULL
      invisible(self)
    },
    
    sim_stateprobs = function(t){
      if(is.null(self$disprog_)){
        stop("You must first simulate disease progression using '$sim_disease'.",
            call. = FALSE)
      }
      private$.stateprobs_ <- indiv_ctstm_sim_stateprobs(self$disprog_,
                                                         self$trans_model,
                                                         t = t)
      invisible(self)
    },
    
    sim_qalys = function(dr = .03, type = c("predict", "random"),
                         lys = TRUE,
                         by_patient = FALSE){
      if(!inherits(self$utility_model, "StateVals")){
        stop("'utility_model' must be an object of class 'StateVals'",
          call. = FALSE)
      }
      type <- match.arg(type)
      private$.qalys_ <- private$sim_wlos(list(self$utility_model), dr, "qalys", type, by_patient,
                                          lys = lys)
      invisible(self)
    },
    
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
      private$.costs_ <- private$sim_wlos(self$cost_models, dr, "costs", type, by_patient, max_t)
      invisible(self)
    },
    
    summarize = function() {
      if (is.null(self$costs_)) {
        stop("Cannot summarize costs without first simulating 'costs_' with '$sim_costs()'.",
               call. = FALSE)
      }
      
      if (is.null(self$qalys_)) {
        stop("Cannot summarize QALYs without first simulating 'qalys_' with '$sim_qalys()'.",
              call. = FALSE)
      }      
      
      stat <- mean
      
      # Costs
      if ("patient_id" %in% colnames(self$costs_)){
        costs_summary <- self$costs_[, lapply(.SD, stat), by = c("category", "dr", "sample", "strategy_id", "state_id"),
                                       .SDcols = "costs"]  
        costs_summary <- costs_summary[, lapply(.SD, sum), by = c("category", "dr", "sample", "strategy_id"),
                                      .SDcols = "costs"]        
      } else{
        costs_summary <- self$costs_[, lapply(.SD, sum), by = c("category", "dr", "sample", "strategy_id"),
                                        .SDcols = "costs"]
      }

      costs_total <- costs_summary[, .(costs = sum(costs)), by = c("dr", "sample", "strategy_id")]
      costs_total[, category := "total"]
      costs_summary <- rbind(costs_summary, costs_total)
      
      # QALYs
      if ("patient_id" %in% colnames(self$qalys_)){
        qalys_summary <- self$qalys_[, lapply(.SD, stat), by = c("dr", "sample", "strategy_id", "state_id"),
                                       .SDcols = "qalys"]  
        qalys_summary <- qalys_summary[, lapply(.SD, sum), by = c("dr", "sample", "strategy_id"),
                                   .SDcols = "qalys"]
      } else{
        qalys_summary <- self$qalys_[, lapply(.SD, sum), by = c("dr", "sample", "strategy_id"),
                                   .SDcols = "qalys"]
      }      

      # Combine
      ce <- list(costs = costs_summary, qalys = qalys_summary)
      class(ce) <- "ce"
      return(ce)
    }
    
  ) # end public
) # end class
