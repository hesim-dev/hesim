# Check ------------------------------------------------------------------------
#' Input validation for class objects
#' 
#' \code{check} is a generic function for validating the inputs of class objects.
#' @param object object to check.
#' @param inner_class When checking a list of objects, the class of elements within
#' the inner most list.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return If validation is successful, returns the object in question; otherwise,
#' informs the user that an error has occurred.  
#' @keywords internal
check <- function (object, ...) {
  UseMethod("check")
}

check_is_class <- function(object, class, name = NULL){
  if (is.null(name)){
    name <- class
  }
  if (!inherits(object, class)){
    stop(paste0("'", name, "' must be of class '", class, "'."),
         call. = FALSE)
  }  
}

check_scalar <- function(x, name){
  if(!(is.atomic(x) && length(x) == 1L)){
    stop(paste0(name, " must be a scalar."))
  }
}

check_dr <- function(dr){
  if(any(table(dr) > 1)){
    stop("You cannot specify the same discount rate twice.",
         call. = FALSE)
  }  
}

check.array <- function(coefs){
  # 'coefs' must be a 3D array (and this has been checked in get_n_samples())
  
  # There are currently no other checks
}

check_StateVals <- function(statevalmods, object, 
                            object_name = c("stateprobs", "disprog")) {
  
  object_name <- match.arg(object_name)
  
  if (!is.null(attr(object, "size"))) {
    expected_size <- attr(object, "size")
  } else {
    stop("'size' attribute missing from ", object_name, ".")
  }
  
  for (i in 1:length(statevalmods)){ # Loop over state value models to check size 
    
    check_size <- function(actual, expected, z = NULL, msg = NULL) {
      if (actual != expected) {
        if (is.null(msg)) {
          msg <- paste0("The number of ", z, " in each 'StateVals' model must equal ",
                        "the number of ", z,  " in the '", object_name, "' object, ",
                        "which is ", expected, ".")
        }
        stop(msg, call. = FALSE)
      }
    }
    
    check_size(get_n_samples(statevalmods[[i]]$params),
               expected_size[["n_samples"]],
               z = "samples")
    check_size(get_id_object(statevalmods[[i]])$n_strategies, 
               expected_size[["n_strategies"]],
               z = "strategies")
    check_size(get_id_object(statevalmods[[i]])$n_patients, 
               expected_size[["n_patients"]],
               z = "patients")
    check_size(
      get_id_object(statevalmods[[i]])$n_states + 1, 
      expected_size[["n_states"]],
      msg = paste0("The number of states in each 'StateVals' model ", 
                   "must be one less (since state values are not applied to the ",
                   "death state) than the number of states in the ",
                   "'", object_name, "' object",
                   ", which is ", expected_size[["n_states"]], ".")
    )
  } # End loop over state value models
}

# Additional utility methods ---------------------------------------------------
#' Absorbing states
#' 
#' Returns a vector of absorbing states from a transition matrix.
#' @param trans_mat A transition matrix in the format from the [`mstate`][mstate::mstate] package. 
#' See [`IndivCtstmTrans`].
#' @keywords internal
absorbing <- function(trans_mat){
  which(apply(trans_mat, 1, function(x) all(is.na(x))))
}

#' Form a list from \code{...}
#' 
#' Form a list of objects from \code{...}.
#' @param ... Objects used to form a list.
#' @return A list of objects from \code{...}.
#' @keywords internal
create_object_list <- function(...){
  objects <- list(...)
  if(length(objects) == 1 & inherits(objects[[1]], "list")){
    objects <- objects[[1]]
  }
  return(objects)
}

# Create list of objects
check_object_list <- function(x, inner_class){
  for (i in 1:length(x)){
    if(!inherits(x[[i]], inner_class)){
      msg <- paste0("Each element in list must be of class '", inner_class, "'")
      stop(msg, call. = FALSE)
    }
  } 
  return(x)
}

new_object_list <- function(..., new_class){
  objects <- create_object_list(...)
  class(objects) <- new_class
  return(objects)
}

object_list <- function(..., inner_class, new_class){
  res <- new_object_list(..., new_class = new_class)
  check_object_list(res, inner_class)
}

# Join objects at specified time points
check_joined_object <- function(x, inner_class, model_list){
  check_object_list(x$models, inner_class)
  
  if(model_list == FALSE){
     check_joined_times(x$models, x$times)
  } else {
    if(!is.list(x$times)){
      stop("'times' must be a list.", call. = FALSE)
    }
    for (i in 1:length(x$times)){
      check_joined_times(x$models[[i]], x$times[[i]])
    }
  } 
  return(x)
}

new_joined_object <- function(..., times, new_class){
  objects <- create_object_list(...)
  res <- list(models = objects, times = times)
  class(res) <- new_class
  return(res)
}


joined_object <- function(..., times, inner_class, new_class, model_list = FALSE){
  res <- new_joined_object(..., times = times, new_class = new_class)
  check_joined_object(res, inner_class, model_list)
}

check_joined_times <- function(objects, times){
  stopifnot(is.vector(times))
  stopifnot(is.numeric(times))
  stopifnot(!is.unsorted(times))
  if(length(objects) != (length(times) + 1)){
    stop("Length of joined models must equal 'times' + 1.",
         call. = FALSE)
  }
}

# list to array
list_to_array <- function(L){
  if (is.matrix(L[[1]]) == TRUE){
      array(unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)))
  } else if (is.vector(L[[1]]) == TRUE){
      array(unlist(L), dim = c(1, length(L[[1]]), length(L)))
  } else{
      stop("List must contain matrices or vectors")
  }
}

# List depth
list_depth <- function(list) {
  ifelse(is.list(list), 1L + max(sapply(list, list_depth)), 0L)
}

# Flatten a nested list
flatten_lists <- function(x) {
  if (!inherits(x, "list")) return(list(x))
  else return(unlist(c(lapply(x, flatten_lists)), recursive = FALSE))
}

# Sample from a posterior distribution
sample_from_posterior <- function(n, n_samples){
  if (n < n_samples){
    return(sample.int(n_samples, n, replace = FALSE))
  } else if (n > n_samples) {
    warning(paste0("The number of requested draws for the probabilistic ",
                   "sensitivity analysis (PSA), 'n', is larger than the number ",
                   "of previously sampled values from the probability ",
                   "distribution of interest. Samples for the PSA have ",
                   "consequently been drawn with replacement."),
            call. = FALSE)
    return(sample.int(n_samples, n, replace = TRUE))
  } else{
    return(1:n)
  }
}

is_whole_number <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

is_1d_vector <- function(x){
  is.atomic(x) && length(dim(x)) <= 1
}

is_3d_array <- function(x){
  is.atomic(x) && length(dim(x)) == 3
}

check_patient_wt <- function(object, result){
  if (is.null(get_id_object(object)[["patient_wt"]])) {
    result[, ("patient_wt") := NULL]
  }
}

format_costs <- function(x, digits){
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

format_qalys <- function(x, digits){
  formatC(x, format = "f", digits = digits)
}

ci_alpha <- function(prob) {
  if (prob > 1 | prob < 0){
    stop("'prob' must be in the interval (0,1)",
         call. = FALSE)
  }
  lower <- (1 - prob)/2
  upper <- 1 - lower
  return(list(lower = lower, upper = upper))
}

format_ci <- function(est, lower, upper, costs = TRUE, digits){
  if (costs){
    est <- format_costs(est, digits = digits)
    lower <- format_costs(lower, digits = digits)
    upper <- format_costs(upper, digits = digits)
  } else{
    est <- format_qalys(est, digits = digits)
    lower <- format_qalys(lower, digits = digits)
    upper <- format_qalys(upper, digits = digits)
  }
  paste0(est, " (",lower, ", ", upper, ")")
}

format_summary_default <- function(x, pivot_from, id_cols, drop_grp) {
  
  if (!is.null(pivot_from)) {
    rhs <- pivot_from
    lhs <- setdiff(id_cols, pivot_from)
    f <- paste(paste(lhs, collapse=" + "), paste(rhs, collapse = " + "),  sep=" ~ ")
    x <- dcast(x, f, value.var = "value", sep = ", ")
  }
  
  # Drop group if desired
  if (drop_grp && ("grp" %in% colnames(x))) {
    n_grps <- length(unique(x$grp))
    if (n_grps == 1) x[, ("grp") := NULL]
  }
  
  # Return
  return(x)
}

# List of matrices -------------------------------------------------------------
matlist <- function(x){
  class(x) <- "matlist"
  return(x)
}

check.matlist <- function(coefs){
  # 'coefs' must be a list (and this has been cheked in get_n_samples())
  
  # Each element of 'coefs' must be a matrix
  matrix_bool <- unlist(lapply(coefs, is.matrix))
  if(sum(!matrix_bool) > 0){
    stop("'coefs' must be a list of matrices.",
         call. = FALSE)
  } 
  
  # Number of rows in each matrix element of 'coefs' must be equal
  coefs_nrows <- unlist(lapply(coefs, nrow))
  if(!all(coefs_nrows[[1]] == coefs_nrows)){
    stop("Number of rows in all 'coefs' matrices must be equal.",
         call. = FALSE)
  } 
}

# Get number of samples --------------------------------------------------------
get_n_samples <- function (x) {
  UseMethod("get_n_samples")
}

get_n_samples.default <- function(x) {
  if ("n_samples" %in% names(x)) {
    return(x[["n_samples"]])
  } else if (is.list(x) && ("n_samples" %in% names(x[[1]]))) {
    return (x[[1]][["n_samples"]])
  } else {
    stop("Could not find 'n_samples'.")
  }
} 

get_n_samples.matlist <- function(x){
  stopifnot(is.list(x))
  if (!is.matrix(x[[1]])){
    stop("'coefs' must be a list of matrices.", call. = FALSE)
  }
  return(nrow(x[[1]]))
}

get_n_samples.array <- function(x){
  stopifnot(is_3d_array(x))
  return(dim(x)[1])
}
