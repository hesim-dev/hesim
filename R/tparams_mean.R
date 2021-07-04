#' Predicted means
#' 
#' Create a list containing means predicted from a statistical model.
#' @param value  Matrix of samples from the distribution of the 
#' mean. Columns denote random samples and rows denote means for different observations.
#' @param ... Arguments to pass to [id_attributes]. Each row in
#' `value` must be a prediction for a `strategy_id`,
#'  `patient_id`, `state_id`, and optionally `time_id` combination.
#'  
#' @return An object of class `tparams_mean`, which is a list containing `value`,
#' `n_samples`, and the ID attributes passed to [id_attributes].
#'  
#' @seealso [tparams]
#' @examples 
#' tparams_mean(value = matrix(1:8, nrow = 4),
#'              strategy_id = rep(1:2, each = 2),
#'              n_strategies = 2,
#'              patient_id = rep(1, 4),
#'              n_patients = 1,
#'              state_id = rep(1:2, times = 2),
#'              n_states = 2)
#'
#' @export
tparams_mean <- function(value, ...){
  stopifnot(is.matrix(value))
  check(new_tparams_mean(value, n_samples = ncol(value), ...),
        ...)
}

new_tparams_mean <- function(value, n_samples, ...){
  l <- c(list(value = value,
              n_samples =  n_samples),
         do.call("new_id_attributes", list(...)))
  class(l) <- "tparams_mean"
  return(l)
}

#' @rdname check
check.tparams_mean <- function(object, ...){
  id_args <- list(...)
  check(do.call("new_id_attributes", id_args))
  for (v in c("strategy_id", "patient_id", "state_id")){
    if (nrow(object$value) != length(id_args[[v]])){
      stop("The length of each ID variable must equal the number of rows in 'value'.",
           call. = FALSE)
    }
  }
  return(object)  
}