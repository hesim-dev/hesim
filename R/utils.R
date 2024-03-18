# Check ------------------------------------------------------------------------
#' Input validation for class objects
#' 
#' \code{check} is a generic function for validating the inputs of class objects.
#' @param object object to check.
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
  if(anyDuplicated(dr) != 0L){
    stop("You cannot specify the same discount rate twice.",
         call. = FALSE)
  }  
}

check_StateVals <- function(models, object, 
                            object_name = c("stateprobs", "disprog")) {
  
  if(!is.list(models)){
    stop("'models' must be a list", call. = FALSE)
  }
  object_name <- match.arg(object_name)
  
  # Generic check for state value model
  for (i in 1:length(models)){
    models[[i]]$check()
  }
  
  # Check that state value models have correct size
  if (!is.null(attr(object, "size"))) {
    expected_size <- attr(object, "size")
  } else {
    stop("'size' attribute missing from ", object_name, ".")
  }
  
  ## The expected number of states depends on the number of absorbing states
  expected_states <- expected_size[["n_states"]]
  if (!is.null(attr(object, "absorbing"))) {
    expected_states <- expected_states - length(attr(object, "absorbing"))
  }
  
  ## Loop over models
  for (i in 1:length(models)){ 
    
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
    
    check_size(get_n_samples(models[[i]]$params),
               expected_size[["n_samples"]],
               z = "samples")
    check_size(get_id_object(models[[i]])$n_strategies, 
               expected_size[["n_strategies"]],
               z = "strategies")
    check_size(get_id_object(models[[i]])$n_patients, 
               expected_size[["n_patients"]],
               z = "patients")
    check_size(
      get_id_object(models[[i]])$n_states, 
      expected_states,
      msg = paste0("The number of states in each 'StateVals' model ", 
                   "must equal the number of states in the ",
                   "'", object_name, "' object less the number of ",
                   "absorbing states, which is ", expected_states, ".")
    )
  } # End loop over state value models
}

# Absorbing states -------------------------------------------------------------
#' Absorbing states
#' 
#' This is a generic function that returns a vector of absorbing states.
#' @param x An object of the appropriate class. When `x` is a `matrix`, 
#' it must be a transition matrix in the format from the 
#' [`mstate`][mstate::mstate] package (see also [`IndivCtstmTrans`]). 
#' @keywords internal
absorbing <- function(x, ...) {
  UseMethod("absorbing")
}

#' @rdname absorbing
absorbing.matrix <- function(x, ...){
  which(apply(x, 1, function(z) all(is.na(z))))
}

#' @rdname absorbing
absorbing.tparams_transprobs <- function(x, ...){
  n_states <- ncol(x$value)
  m <- apply(x$value, c(1, 2), mean)
  sum_zero <- apply(m, 1, function (z) sum(z == 0))
  absorbing <- which(sum_zero == (n_states - 1))
  if (length(absorbing) == 0) absorbing <- NULL
  absorbing
}

# Additional utility methods ---------------------------------------------------
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

default_colnames <- function(x) {
  paste0("x", 1:ncol(x))
}

# Summarize parameter values ---------------------------------------------------
summarize_params <- function(x, ...) {
  UseMethod("summarize_params")
}

summarize_params.numeric <- function(x, probs, param_name = "param",
                                     param_values, ...) {
  stats <- as.data.table(t(c(
    mean = mean(x),
    sd = stats::sd(x),
    stats::quantile(x, probs = probs)
  )))
  y <- data.table(
    param = param_values,
    stats
  )
  setnames(y, "param", param_name)
  y
}

summarize_params.integer <- function(x, probs, param_name = "param",
                                     param_values, ...) {
  summarize_params.numeric(x, probs, param_name, param_values)
}

summarize_params.matrix <- function(x, probs, param_name = "param",
                                    param_values = NULL,...) {
  
  apply_quantile <- function(x, probs) {
    y <- apply(x, 2, stats::quantile, probs = probs)
    if (length(probs) == 1) {
      y <- matrix(y, ncol = 1)
      colnames(y) <- paste0(probs * 100, "%")
      return(y)
    } else{
      return(t(y))
    }
  }

  if (is.null(param_values)) param_values <- colnames(x)
  x <- data.table(
    param = param_values,
    mean = apply(x, 2, mean),
    sd = apply(x, 2, stats::sd),
    apply_quantile(x, probs)
  )
  setnames(x, "param", param_name)
  x
}

summarize_params.data.table <- function(x, probs, param_name = "param",
                                        param_values = NULL, cols = NULL,...) {
  
  if (is.null(cols)) cols <- colnames(x)
  x <- as.matrix(x)
  summarize_params.matrix(x, probs = probs, param_name = param_name,
                          param_values = param_values)
}  

# List of matrices -------------------------------------------------------------
coeflist <- function(coefs){
  if (inherits(coefs, "list")) {
    coefs <- lapply(coefs, function(x) {
      x <- as.matrix(x)
      if (is.null(colnames(x))) colnames(x) <- default_colnames(x)
      x
    })
  } else {
    stop("'coefs' must be a list.", call. = FALSE)
  }
  class(coefs) <- "coeflist"
  return(coefs)
}

check.coeflist <- function(object, ...){
  coefs = object
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

get_n_samples.coeflist <- function(x){
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

#' Code to use the hesim package inline. Not directly called by the user.
#' @param ... arguments
#' @examples
#' library(Rcpp)
#' sourceCpp(code="
#' // [[Rcpp::depends(hesim)]]
#' // [[Rcpp::depends(RcppArmadillo)]]
#' #include <hesim/stats/distributions.h>
#' // [[Rcpp::export]]
#' double test_inline_gengamma(double mu, double sigma, double Q) {
#'   hesim::stats::gengamma gg(mu, sigma, Q);
#'   return gg.random();
#' }")
#' set.seed(12345)
#' test_inline_gengamma(1.0, 1.0, 1.0)

#' @rdname plugin
inlineCxxPlugin <- function(...) {
    ismacos <- Sys.info()[["sysname"]] == "Darwin"
    openmpflag <- if (ismacos) "" else "$(SHLIB_OPENMP_CFLAGS)"
    plugin <- Rcpp::Rcpp.plugin.maker(include.before = "#include <hesim.h>",
                                      libs           = paste(openmpflag,
                                                             "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)"),
                                      package        = "hesim")
    settings <- plugin()
    settings$env$PKG_CPPFLAGS <- paste("-I../inst/include", openmpflag)
    ## if (!ismacos) settings$env$USE_CXX11 <- "yes"
    settings
}
