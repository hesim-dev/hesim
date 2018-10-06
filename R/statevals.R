# stateval_means ----------------------------------------------------------------
#' Mean values for health states
#' 
#' Estimated means that are used to simulate values assigned to health states with
#'  \code{\link{StateVals}}.
#' @param values Typically a matrix where each column denotes a health state and each
#' row denotes a random sample of the value assigned to that health state for the
#' probabilistic sensitivitiy analysis. If health state values vary by strategy, then
#' it can be an array of matrices where each matrix denotes the values associated with
#' a different treatment strategy. 
#' @param strategy_id If a matrix, a vector denoting the strategy ID's to be modeled; if an array,
#' a vector denoting the strategy ID associated with each matrix in the array.
#' @param patient_id The patient ID's that will be modeled in the simulation. Note that state
#' values are assumed to be constant across patients when using objects of class "stateval_means". 
#' 
#' @return An object of class "stateval_means", which is a list containing \code{values},
#' \code{strategy_id}, and \code{patient_id}.
#' @examples 
#' # Cost estimates in 2 health states for a model with 2 treatment strategies and 3 patients
#' gamma_params <- mom_gamma(c(5000, 7000), c(1000, 1200))
#' n <- 3
#' vals <- matrix(rgamma(2 * n, 
#'                       shape = gamma_params$shape, 
#'                       scale = gamma_params$scale),
#'                nrow = n, ncol = 2, byrow = TRUE)
#' stval_ests <- stateval_means(values = vals,
#'                            strategy_id = c(1, 2),
#'                            patient_id = c(1, 2, 3))
#' print(stval_ests)
#' stateval_mod <- create_StateVals(stval_ests)
#' head(stateval_mod$sim(t = c(1, 2, 3), type = "predict"))
#'
#' @export
stateval_means <- function(values, strategy_id, patient_id){
  if(!(is.matrix(values) | is.array(values))){
    stop("'values' must be a matrix or an array.", 
         call. = FALSE)
  }
  if (!is.matrix(values)){
    if (length(strategy_id) != dim(values)[3]){
      stop("The length of 'strategy_id' must equal the number of matrices in 'values'",
           call. = FALSE)
    }
  }
  l <- list(values = values, strategy_id = strategy_id,
            patient_id = patient_id)
  class(l) <- "stateval_means"
  return(l)
}

# StateVals --------------------------------------------------------------------
#' Create \code{StateVals} object
#' 
#' \code{create_StateVals} is a generic function for creating an object of class
#'  \code{\link{StateVals}} from a fitted statistical model or \code{\link{stateval_means}}
#'  object. 
#' @param object A model object of the appropriate class.
#' @param data An object of class "expanded_hesim_data" returned by 
#' \code{\link{expand.hesim_data}}. Must be expanded by the data tables "strategies",
#' "patients", and "states".
#' @param n Number of random observations of the parameters to draw when parameters 
#' are fit using a statistical model.
#' @param point_estimate If \code{TRUE}, then the point estimates are returned and and no samples are drawn.
#' @param ... Further arguments passed to or from other methods. Currently unused. 
#' @return Returns an \code{\link{R6Class}} object of class \code{\link{StateVals}}.
#' @seealso \code{\link{StateVals}}
#' @export
create_StateVals <- function(object, ...){
  UseMethod("create_StateVals", object)
} 
 
#' @rdname create_StateVals
#' @export  
create_StateVals.lm <- function(object, data = NULL, n = 1000,
                                point_estimate = FALSE, ...){
  params <- create_params(object, n, point_estimate) 
  input_mats <- create_input_mats(object, data)
  return(StateVals$new(input_mats = input_mats, params = params))
}

#' @rdname create_StateVals 
#' @export
create_StateVals.stateval_means <- function(object, ...){
  params <- create_params(object)
  input_mats <- create_input_mats(object)
  return(StateVals$new(input_mats = input_mats, params = params))
}

# Manual documentation in StateVals.Rd
#' @export
StateVals <- R6::R6Class("StateVals",
  public = list(
    input_mats = NULL,
    params = NULL,

    initialize = function(input_mats, params) {
      self$input_mats <- input_mats
      self$params <- params
    },
    
    sim = function(t, type = c("predict", "random")){
      type <- match.arg(type)
      self$check()
      res <- data.table(C_statevals_sim(self, t, type))
      res[, sample := sample + 1]
      return(res[])
    },
    
    check = function(){
      if(!inherits(self$input_mats, "input_mats")){
        stop("'data' must be an object of class 'input_mats'",
            call. = FALSE)
      }
      if(!inherits(self$params, c("params_mean", "params_lm"))){
          stop("Class of 'params' is not supported. See documentation.",
               call. = FALSE)
      }
    }
  )
)

