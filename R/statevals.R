# PartSurvStateVals ------------------------------------------------------------
#' Form \code{StateVals} object
#' 
#' \code{create_StateVals} is a generic function for forming an object of class
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

