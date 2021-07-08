# tparams_transprobs -----------------------------------------------------------
#' Transition probabilities
#' 
#' Create a list containing predicted transition probabilities at discrete times.
#' Since the transition probabilities have presumably already been predicted 
#' based on covariate values, no input data is required for
#' simulation. The class can be instantiated from either an `array`, 
#' a `data.table`, a `data.frame`, or a [`tpmatrix`]. This is the object in
#' `hesim` used to specify the transition probabilities required to simulate 
#' Markov chains with the [`CohortDtstmTrans`] class. 
#' 
#' @param object An object of the appropriate class. 
#' @param tpmatrix_id An object of class [`tpmatrix_id`] (or an equivalent 
#' `data.table` with the same ID columns as returned by `tpmatrix_id()`).
#' @param ... Further arguments passed to or from other methods. Currently unused. 
#' @param times An optional numeric vector of distinct times to pass to
#'  [time_intervals] representing time intervals indexed by the 4th dimension of 
#'  the array. May either be the start or the end of intervals. 
#'  This argument is not required if there is only one time interval.
#' @param grp_id An optional numeric vector of integers denoting the subgroups. Must
#' be the same length as the 3rd dimension of the array.
#' @param patient_wt An optional numeric vector denoting the weight to apply to each 
#' patient within a subgroup. Must be the same length as the 3rd dimension of the array.
#' 
#' @details The format of `object` depends on its class: 
#' \describe{
#' \item{array}{
#'  Either a 3D or a 6D array is possible.
#'\itemize{
#'  \item If a 3D array, then each slice is a 
#'  square transition probability matrix. In this case
#' `tpmatrix_id` is required and each matrix slice corresponds to the same 
#'  numbered row in `tpmatrix_id`. The number of matrix slices must equal the number 
#'  of rows in `tpmatrix_id`.
#'  
#'  \item If a 6D array, then the dimensions of the array should be indexed as follows: 
#'  1st (`sample`), 2nd (`strategy_id`), 3rd (`patient_id`),
#'  4th (`time_id`), 5th (rows of transition matrix), and
#'  6th (columns of transition matrix). In other words, an index of
#'  `[s, k, i, t]` represents the transition matrix for the `s`th 
#'  sample, `k`th treatment strategy, `i`th patient, and `t`th
#'  time interval.
#' }
#' }
#'  
#'  \item{data.table}{Must contain the following:
#'  \itemize{
#'  \item ID columns for the parameter sample (`sample`), 
#'  treatment strategy (`strategy_id`), and patient (`patient_id`).
#'  If the number of time intervals is greater than 1 it must also contain the
#'  column `time_start` denoting the starting time of a time interval. A column 
#'  `patient_wt` may also be used to denote the weight to apply to each
#'  patient.
#'  \item Columns for each element of the transition probability matrix. 
#'  They should be prefixed with "prob_" and ordered rowwise. 
#'  For example, the following columns would be used for a 2x2 transition
#'   probability matrix:
#'  `prob_1` (1st row, 1st column), 
#'  `prob_2` (1st row, 2nd column), 
#'  `prob_3` (2nd row, 1st column), and 
#'  `prob_4` (2nd row, 2nd column).
#'  }
#'  }
#'  
#'  \item{data.frame}{Same as `data.table`.}
#'  
#'  \item{tpmatrix}{An object of class [`tpmatrix`].}
#' }
#' 
#' A `tparams_transprobs` object is also instantiated when creating a
#' cohort discrete time state transition model using [`define_model()`].
#' 
#' 
#' @return An object of class `tparams_transprobs`, 
#' which is a list containing `value` and relevant ID attributes. The element 
#' `value` is an array of predicted transition probability matrices from the 
#' probability distribution of the underlying statistical model. Each matrix in 
#' `value` is a prediction for a `sample`, `strategy_id`, `patient_id`, and 
#' optionally `time_id` combination.
#'
#' @seealso A `tparams_transprobs` object is used to store the "parameters" of 
#' the transition component of a cohort discrete time state transition 
#' model (cDTSTM). You can create such an object with `CohortDtstmTran$new()`.
#' 
#' [tpmatrix()] and [tpmatrix_id()] provide a convenient way to construct a 
#' `tparams_transprobs` object in a flexible way. [`define_model()`] is, in turn,
#' a convenient way to construct a [`tpmatrix`] object using mathematical 
#' expressions; in this case, an entire cDTSTM can be instantiated from a model 
#' definition using [create_CohortDtstm.model_def()]. Detailed examples
#' are provided in `vignette("markov-cohort")` and 
#' `vignette("markov-inhomogeneous-cohort")`
#' 
#' The output of a `tparams_transprobs` object is rather verbose. It can be
#' helpful to check the output by converting it to a `data.table` (containing
#' both the ID variables and flattened transition probability matrices)
#' with [as.data.table.tparams_transprobs()]. Transition probabilities can
#' also be summarized (across parameter samples) using 
#' [summary.tparams_transprobs()].
#' 
#' 
#' @example man-roxygen/example-tparams_transprobs.R
#' @rdname tparams_transprobs
#' @export
tparams_transprobs <- function(object, ...){
  if (missing(object)){
    stop("'object' is missing with no default.")
  }
  UseMethod("tparams_transprobs", object)
} 

new_tparams_transprobs <- function(object, ...){
  UseMethod("new_tparams_transprobs", object)
} 

create_tparams_transprobs <- function(value, ...){
  C_normalize_transprobs(value)
  l <- c(list(value = value),
         do.call("new_id_attributes", list(...)))
  class(l) <- "tparams_transprobs"
  return(l)
}

#' @rdname check
check.tparams_transprobs <- function(object){
  stopifnot(is.array(object$value))
  stopifnot(is.numeric(object$sample))
  stopifnot(is.numeric(object$n_samples))
  id_args <- object[names(object) != "value"]
  check(do.call("new_id_attributes", id_args))
  return(object)
}

tparams_transprobs_id <- function(x) {
  # ID attributes
  id_args <- list()
  id_names <- c("sample", "strategy_id", "patient_id")
  size_names <- c("n_samples", "n_strategies", "n_patients")
  for (i in 1:length(id_names)){
    id_args[[id_names[i]]] <- x[[id_names[i]]]
    id_args[[size_names[i]]] <- length(unique(x[[id_names[i]]]))
  }
  
  ## Time interval
  if (!is.null(x[["time_start"]])){
    time_intervals <- time_intervals(unique(x[["time_start"]])) 
    pos <- match(x[["time_start"]], time_intervals$time_start)
    id_args[["time_id"]] <- time_intervals$time_id[pos]
    id_args[["time_intervals"]] <- time_intervals
    id_args[["n_times"]] <- length(time_intervals$time_id)
  } else{
    id_args[["time_id"]] <- rep(1, length(id_args$sample))
    id_args[["time_intervals"]] <- data.table(time_id = 1,
                                              time_start = 0,
                                              time_stop = Inf)
    id_args[["n_times"]] <- 1
  }
  
  # Group ID and patient weight
  id_args[["grp_id"]] <- x[["grp_id"]]
  id_args[["patient_wt"]] <- x[["patient_wt"]]
  return(id_args)
}

new_tparams_transprobs_array6 <- function (object, times = NULL, 
                                           grp_id = NULL, patient_wt = NULL) {
  
  # Reshape array
  dims <- c(dim(object)[5], dim(object)[6], prod(dim(object)[1:4]))
  value <- array(c(aperm(object, perm = c(5, 6, 4, 3, 2, 1))),
                 dim = dims)
  
  # ID attributes
  n_patients <- 1:dim(object)[3]
  id_df <- expand.grid(time_id = 1:dim(object)[4],
                       patient_id = n_patients,
                       strategy_id = 1:dim(object)[2],
                       sample = 1:dim(object)[1])
  if (!is.null(grp_id) | !is.null(patient_wt)){
    check_var <- function(x, name){
      if (!is.null(x) & !length(x) %in% c(1, n_patients)){
        stop(paste0("The length of '", name, "' must be equal to the 3rd dimension ",
                    "of the array (i.e., the number of patients).", call. = FALSE),
             call. = FALSE)
      }
    }
    check_var(grp_id, "grp_id")
    check_var(patient_wt, "patient_wt")
    tmp_dt <- data.table(patient_id = 1:dim(object)[3],
                         grp_id = grp_id,
                         patient_wt = patient_wt)
    id_df <- cbind(id_df, 
                   tmp_dt[match(id_df$patient_id, tmp_dt$patient_id),
                          list(grp_id, patient_wt)])
  }
  n_df <- list(n_samples = dim(object)[1],
               n_strategies = dim(object)[2],
               n_patients = dim(object)[3],
               n_times = dim(object)[4])
  id_args <- list()
  for (v in colnames(id_df)){
    id_args[[v]] <- id_df[[v]]
  }
  for (v in names(n_df)){
    id_args[[v]] <- n_df[[v]]
  }
  if (n_df$n_times == 1 ){
    id_args$time_intervals <- time_intervals(0)
  } else{
    if (is.null(times)){
      stop(paste0("'times' cannot be NULL if the number of time ",
                  "intervals is greater than 1"), call. = FALSE)
    }
    id_args$time_intervals <- time_intervals(times)
  }
  
  # Return
  return(list(value = value, id_args = id_args))
}


new_tparams_transprobs.array <- function (object, tpmatrix_id = NULL, times = NULL, 
                                          grp_id = NULL, patient_wt = NULL) {
  # Checks
  n_dim <- length(dim(object))
  if(!n_dim %in% c(3, 6)){
    stop("'object' must be a 3D or 6D array.", call. = FALSE)
  }  
  
  # Return
  if (n_dim == 6) {
    y <- new_tparams_transprobs_array6(object, times, grp_id, patient_wt) 
    return(do.call("create_tparams_transprobs", c(list(value = y$value), y$id_args)))
  } else{
    check_is_class(tpmatrix_id, "data.frame", "tpmatrix_id")
    if (nrow(tpmatrix_id) != dim(object)[3]) {
      stop(paste0("The third dimension of the array 'object' must equal the ", 
                  "number or rows in 'tpmatrix_id'."), call. = FALSE)
    }
    id_args <- tparams_transprobs_id(tpmatrix_id)
    return(do.call("create_tparams_transprobs", c(list(value = object), id_args)))
  }
}

#' @rdname tparams_transprobs
#' @export
tparams_transprobs.array <- function (object, tpmatrix_id = NULL, times = NULL, 
                                      grp_id = NULL, patient_wt = NULL, ...) {
  res <- new_tparams_transprobs.array(object = object, tpmatrix_id = tpmatrix_id,
                                      times = times, grp_id = grp_id, 
                                      patient_wt = patient_wt)
  return(check(res))
}

new_tparams_transprobs.data.table <- function (object) {
  id_args <- tparams_transprobs_id(object)
  indices <- grep("^prob_", colnames(object))
  if (length(indices) == 0) {
    stop("No columns with names starting with 'prob_'.")
  }
  prob_mat <- as.matrix(object[, colnames(object)[indices], with = FALSE])
  value <- as_array3(prob_mat)
  return(do.call("create_tparams_transprobs", c(list(value = value), id_args)))
}

#' @rdname tparams_transprobs
#' @export
tparams_transprobs.data.table <- function (object, ...) {
  return(check(new_tparams_transprobs(object)))
}

#' @rdname tparams_transprobs
#' @export
tparams_transprobs.data.frame <- function (object, ...) {
  res <- new_tparams_transprobs(data.table(object))
  return(check(res))
}

#' @rdname tparams_transprobs
#' @export
tparams_transprobs.tpmatrix <- function(object, tpmatrix_id, ...) {
  check_is_class(tpmatrix_id, "data.frame", "tpmatrix_id")
  if (nrow(object) != nrow(tpmatrix_id)) {
    stop("'object' and 'tpmatrix_id' must have the same number of rows.",
         call. = FALSE)
  }
  value <- as_array3(object)
  id_args <- tparams_transprobs_id(tpmatrix_id)
  return(do.call("create_tparams_transprobs", c(list(value = value), id_args)))
}

tparams_transprobs.eval_model <- function(object, ...){
  id_index <- attr(object$tpmatrix, "id_index")
  return(tparams_transprobs(object$tpmatrix, object$id[[id_index]]))
}

# summary.tparams_transprobs ---------------------------------------------------
#' Summarize `tparams_transprobs` object
#' 
#' The `summary()` method summarizes a [`tparams_transprobs`] object containing 
#' predicted transition probabilities; summary statistics are computed for each 
#' possible transition by the relevant ID variables. The `print()` method
#' summarizes the object using `summary.tparams_transprobs()` and prints it to the
#' console. 
#' 
#' @inheritParams summary.tpmatrix
#' @param object,x A [`tparams_transprobs`] object.
#' 
#' @return `summary.tparams_transprobs()` uses the [summary.tpmatrix()]
#' method under the hood and returns the same output. 
#' 
#' @seealso See [`tparams_transprobs`] for an example use of the summary and
#' print methods. Also see [summary.tpmatrix()] for more information on what is 
#' returned.
#' 
#' @export
summary.tparams_transprobs <- function(object, probs = c(0.025, 0.975), 
                                       unflatten = FALSE, ...) {
  # Convert value element to tpmatrix
  object_dt <- as_tbl2(object$value)
  tpmat <- tpmatrix(object_dt, states = colnames(object$value))
  
  # Create tpmatrix_id object
  id_dt <- make_id_data_table(object)
  setattr(id_dt, "class", c("tpmatrix_id", "data.table", "data.frame"))
  
  # Summarize
  summary(tpmat, id = id_dt, probs = probs, unflatten = unflatten)
}

# print.tparams_transprobs -----------------------------------------------------
#' @rdname summary.tparams_transprobs
#' @param ... For the `print()` method, arguments to pass to `summary.tparams_transprobs`;
#' currently unused for the `summary()` method.
#' @export
print.tparams_transprobs <- function(x, ...) {
  cat("A \"tparams_transprobs\" object \n\n")
  cat("Summary of transition probabilities:\n")
  print(summary(x, ...))
  invisible(x)
}

# as.data.table.tparams_transprobs ---------------------------------------------
#' Coerce to `data.table`
#' 
#' Creates a `data.table` that combines the transition probability matrices 
#' and ID columns from a [tparams_transprobs] object. This is often useful for 
#' debugging. 
#' @param x A [tparams_transprobs] object.
#' @param prefix,sep Arguments passed to [tpmatrix_names()] for naming
#' the transition probability columns. The `states` argument is based on
#' the column names (i.e., names of the second dimension) of the `$value`
#' element of `x`; if `NULL`, then states are named `s1`, ..., `sh` where h is 
#' the number of states.
#' @param ... Currently unused. 
#' 
#' @seealso [tparams_transprobs]
#' @return A `data.table` with one row for each transition probability matrix.
#' 
#' @examples 
#' # Create tparams_transprobs object
#' hesim_dat <- hesim_data(strategies = data.frame(strategy_id = 1:2),
#'                         patients = data.frame(patient_id = 1:3))
#' input_data <- expand(hesim_dat, by = c("strategies", "patients"))    
#' tpmat_id <- tpmatrix_id(input_data, n_samples = 2)      
#' p_12 <- runif(nrow(tpmat_id), .6, .7) + 
#'   .05 * (tpmat_id$strategy_id == 2)
#' tpmat <- tpmatrix(
#'   C, p_12,
#'   0, 1
#' )
#' tprobs <- tparams_transprobs(tpmat, tpmat_id)
#' 
#' # Convert to data.table
#' as.data.table(tprobs)
#' as.data.table(tprobs, prefix = "")
#' as.data.table(tprobs, prefix = "", sep = ".")
#' 
#' @export
as.data.table.tparams_transprobs <- function(x, ..., prefix = "prob_", sep = "_"){
  probs <- as_tbl2(x$value, prefix = prefix, sep = sep)
  id_dt <- as.data.table(x[c("sample", "strategy_id", "patient_id")])
  time_dt <- x$time_intervals[match(x$time_id, x$time_intervals$time_id)]
  x_dt <- data.table(id_dt, time_dt, probs)
  for (v in c("n_samples", "n_strategies", "n_patients", "n_times")){
    setattr(x_dt, v, x[[v]])
  }
  return(x_dt)
}