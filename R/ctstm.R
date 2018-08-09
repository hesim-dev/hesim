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
#' @param ... Further arguments passed to or from other methods. Currently unused. 
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
  return(CtstmTrans$new(data = input_data, params = params, trans_mat = trans_mat))
}

#' @export
#' @rdname create_CtstmTrans
create_CtstmTrans.flexsurvreg <- function(object, data, trans_mat, n = 1000, point_estimate = FALSE, ...){
  input_data <- create_input_data(object, data, id_vars = c("strategy_id", "patient_id", "transition_id"))
  params <- create_params(object, n = n, point_estimate = point_estimate)
  return(CtstmTrans$new(data = input_data, params = params, trans_mat = trans_mat))
}

CtstmTrans <- R6::R6Class("CtstmTrans",
  private = list(
      summary = function(t, type = c("hazard", "cumhazard")){
        #self$check()
        type <- match.arg(type)
        res <- data.table(C_ctstm_summary(self, t, type))
        res[, trans := trans + 1]
        res[, sample := sample + 1]
        if (type == "hazard") setnames(res, "value", "hazard")
        if (type == "cumhazard") setnames(res, "value", "cumhazard")
        return(res[])
      }
    ), 
                          
  public = list(
    data = NULL,
    params = NULL,
    trans_mat = NULL,
    initialize = function(data, params, trans_mat) {
      self$data <- data
      self$params <- params
      self$trans_mat <- trans_mat
    },
    
    hazard = function(t){
      private$summary(t, "hazard")
    },
    
    cumhazard = function(t){
       private$summary(t, "cumhazard")
    },
    
    pmatrix = function(){
      print("Code goes here.")
    }
  )
)

# CtStmIps ---------------------------------------------------------------------
IndivCtstm <- R6::R6Class("IndivCtstm",
  public = list(
    trans_models = NULL,
    trans_mat = NULL,
    start = NULL,
    max_t = NULL,
    max_age = NULL,
    initialize = function(trans_models = NULL, trans_mat, start = 1, max_t = 100, max_age = NULL) {
      self$trans_models <- trans_models
      self$trans_mat <- trans_mat
      self$start <- start
      self$max_t = max_t
      self$max_age = mage_age
    }
  )
)
