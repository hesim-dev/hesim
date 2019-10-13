# Transition probability matrix class ------------------------------------------
#' @export
`[.transprob_matrix` <- function(x, i, j, ...) {
  n_rows <- n_cols <- sqrt(length(x))
  index <- (i - 1) * n_cols + j
  y <- as.vector(x)
  return(y[index])
}
#' @examples 
#' p <- 1:4
#' attr(p, "class") <- "transprob_matrix"
#' p[1, 1]
#' p[1, 2]
#' p[2, 1]
#' p[2, 2]

# CohortDtstmTrans -------------------------------------------------------------
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
    params = NULL,
    input_mats = NULL,
    start_stateprobs = NULL,
    cycle_length = NULL,
    
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
    
    sim_stateprobs = function(n_cycles){
      times <- seq(0, n_cycles/self$cycle_length, length.out = n_cycles + 1)
      stprobs <- C_cohort_dtstm_sim_stateprobs(self, 
                                               times,
                                               self$params$n_samples)
      stprobs <- data.table(stprobs)
      stprobs[, sample := sample + 1]
      stprobs[, state_id := state_id + 1]
      return(stprobs[])
    }
  )
)

# create_CohortDtstmTrans ------------------------------------------------------
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
#' @export
CohortDtstm <- R6::R6Class("CohortDtstm",
  public = list(
    trans_model = NULL,
    utility_model = NULL,
    cost_models = NULL,
    stateprobs_ = NULL,
    qalys_ = NULL,
    costs_ = NULL,
    initialize = function(trans_model = NULL, utility_model = NULL, cost_models = NULL) {
      self$trans_model <- trans_model
      self$utility_model = utility_model
      self$cost_models = cost_models
    },
    
    sim_stateprobs = function(n_cycles){
      self$stateprobs_ <- self$trans_model$sim_stateprobs(n_cycles)
      setattr(self$stateprobs_, "class", 
              c("stateprobs", "data.table", "data.frame"))
      invisible(self)
    },
    
    sim_qalys = function(dr = .03,
                         method = c("trapz", "riemann_left", "riemann_right")){
      self$qalys_ <- sim_qalys(self$stateprobs_, self$utility_model, dr, method)
      invisible(self)
    },
    
    sim_costs = function(dr = .03, 
                         method = c("trapz", "riemann_left", "riemann_right")){
      self$costs_ <- sim_costs(self$stateprobs_, self$cost_models, dr, method)
      invisible(self)
    },
    
    summarize = function() {
      check_summarize(self)
      return(summarize_ce(self$costs_, self$qalys_))
    }
  )
)
