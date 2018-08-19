# params_statevals -------------------------------------------------------------
#' Parameters for state values
#' 
#' Parameters that are used to simulate values assigned to health states with
#'  \code{\link{StateVals}}.
#' @param values Typically a matrix where each column denotes a health state and each
#' row denotes a random sample of the value assigned to that health state for the
#' probabilistic sensitivitiy analysis. If health state values vary by strategy, then
#' it can be an array of matrices where each matrix denotes the values associated with
#' a different treatment strategy. 
#' @param timefun An R function that can be used to modify \code{values} as a function
#' of time. 
#' 
#' @return An object of class "params_statevals", which is a list containing \code{values} and
#' \code{timefun}.
#' @examples 
#' # Cost estimates in 2 health states
#' gamma_params <- mom_gamma(c(5000, 7000), c(1000, 1200))
#' n <- 1000
#' vals <- matrix(rgamma(2 * n, 
#'                       shape = gamma_params$shape, 
#'                       scale = gamma_params$scale),
#'                nrow = n, ncol = 2, byrow = TRUE)
#' params_stvals <- params_statevals(values = vals)
#' print(params_stvals)
#'
#' @export
params_statevals <- function(values, timefun = NULL){
  if(!is.matrix(values) | !is.array(values)){
    stop("'values' must be a matrix or an array.")
  }
  l <- list(values = values, timefun = timefun)
  class(l) <- "params_statevals"
  return(l)
}

# PartSurvStateVals ------------------------------------------------------------
#' Create \code{StateVals} object
#' 
#' \code{create_StateVals} is a generic function for creating an object of class
#'  \code{\link{StateVals}} from a fitted statistical model. 
#' @param object A fitted statistical model object of the appropriate class. Supports
#' \code{\link{lm}}.
#' @param data An object of class "expanded_hesim_data" returned by 
#' \code{\link{expand_hesim_data}}. Must be expanded by the data tables "strategies",
#' "patients", and "states".
#' @param n Number of random observations of the parameters to draw.
#' @param point_estimate If \code{TRUE}, then the point estimates are returned and and no samples are drawn.
#' @param ... Further arguments passed to or from other methods. Currently unused. 
#' @return Returns an \code{\link{R6Class}} object of class \code{\link{StateVals}}.
#' @seealso \code{\link{StateVals}}
#' @export
create_StateVals <- function(object, data, n = 1000, point_estimate = FALSE){
  if (!inherits(object, c("lm"))){
    stop("Class of 'object' is not supported. See documentation.",
         call. = FALSE)
  }
  input_data <- create_input_data(object, data, id_vars = c("strategy_id", "patient_id", "state_id")) 
  params <- create_params(object, n, point_estimate)
  return(StateVals$new(data = input_data, params = params))
}

# Manual documentation in StateVals.Rd
#' @export
StateVals <- R6::R6Class("StateVals",
  public = list(
    data = NULL,
    params = NULL,

    initialize = function(data, params) {
      self$data <- data
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
      if(!inherits(self$data, "input_data")){
        stop("'data' must be an object of class 'input_data'",
            call. = FALSE)
      }
      if(!inherits(self$params, c("params_lm"))){
          stop("Class of 'params' is not supported. See documentation.",
               call. = FALSE)
      }
    }
  )
)

