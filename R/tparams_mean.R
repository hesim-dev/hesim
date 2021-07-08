# tparams_mean -----------------------------------------------------------------
#' Predicted means
#' 
#' Create a list containing means predicted from a statistical model.
#' 
#' @param value  Matrix of samples from the distribution of the 
#' mean. Columns denote random samples and rows denote means for different 
#' observations.
#' @param ... Arguments to pass to [id_attributes]. Each row in
#' `value` must be a prediction for a `strategy_id`,
#'  `patient_id`, `state_id`, and optionally `time_id` combination.
#'  
#' @note The `tparams_mean()` constructor would not normally be used by users; instead,
#' a `tparams_mean` object is typically created automatically as part of the 
#' [`StateVals`] class with [create_StateVals()].
#'  
#' @return An object of class `tparams_mean`, which is a list containing `value`,
#' `n_samples`, and the ID attributes passed to [id_attributes].
#'  
#' @seealso A `tparams_mean` object is a type of [transformed parameter][tparams]
#' object and is a supported class type of the `params` field of the [`StateVals`]
#' class. See the documentation for [create_StateVals()] and [stateval_tbl()]
#' for examples of how to create`StateVals` objects. Predicted means can be 
#' summarized across parameter samples using [summary.tparams_mean()].
#' 
#' @example man-roxygen/example-tparams_mean.R
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

# summary.tparams_mean ---------------------------------------------------------
#' Summarize `tparams_mean` object
#' 
#' The `summary()` method summarizes a [`tparams_mean`] object containing 
#' predicted means; summary statistics are computed for each 
#' combination of the ID variables. The `print()` method
#' summarizes the object using `summary.tparams_mean()` and prints it to the
#' console. 
#' 
#' @inheritParams summary.tparams_transprobs
#' @param object,x A [`tparams_mean`] object.
#' @param ... Currently unused.
#' 
#' @return A `data.table` with columns for (i) the ID variables, 
#' (ii) the mean of each parameter across parameter samples (`mean`),
#' (iii) the standard deviation of the parameter samples (`sd`), and
#' (iv) quantiles of the parameter samples corresponding to the `probs` argument.
#' 
#' @seealso See [`tparams_mean`] for an example use of the summary and
#' print methods.
#' 
#' @export
summary.tparams_mean <- function(object, probs = c(0.025, 0.975), ...) {
  q <- apply(object$value, 1, stats::quantile, probs = probs)
  if (is.matrix(q)) {
    q <- t(q)
  } else{
    q <- as.matrix(q)
    colnames(q) <- paste0(probs * 100, "%")
  }
  
  data.table(
    make_id_data_table(object),
    mean = apply(object$value, 1, mean),
    sd = apply(object$value, 1, stats::sd),
    q
  )
}

# print.tparams_mean -----------------------------------------------------------
#' @rdname summary.tparams_mean
#' @export
print.tparams_mean <- function(x, ...) {
  cat("A \"tparams_mean\" object \n\n")
  cat("Summary of means:\n")
  print(summary(x, ...))
  invisible(x)
}