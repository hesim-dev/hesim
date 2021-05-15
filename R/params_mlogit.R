coeflist_to_array <- function(coefs) {
  coefs <- lapply(coefs, as.matrix)
  dims <- sapply(coefs, dim)
  if (!all(dims[1, ] == dims[1, 1])) {
    stop("The number of rows in each element of 'coefs' must be the same.",
         .call = FALSE)
  }
  if (!all(dims[2, ] == dims[2, 1])) {
    stop("The number of columns in each element of 'coefs' must be the same.",
         .call = FALSE)
  }
  
  array(unlist(coefs), 
        dim = c(nrow(coefs[[1]]), ncol(coefs[[1]]), length(coefs)),
        dimnames = list(NULL, colnames(coefs[[1]]), names(coefs)))
}


# params_mlogit() --------------------------------------------------------------
#' Parameters of a multinomial logit model
#' 
#' Store the parameters of a fitted multinomial logistic 
#' regression model. The model is used to predict probabilities of \eqn{K} 
#' classes.
#' @param coefs  A 3D array of stacked matrices. The number of matrices (i.e.,
#' the number of slices in the cube) should be equal to \eqn{K-1}. Each 
#' matrix contains samples of the regression coefficients under sampling uncertainty
#' corresponding to a particular class. Rows index parameter samples and 
#' columns index coefficients.
#' 
#' @details Multinomial logit models are used to predict the probability of 
#' membership for subject \eqn{i} in each of \eqn{K} classes as a function of covariates:
#' \deqn{Pr(y_i = c) = \frac{e^{\beta_c x_i}}{\sum_{k=1}^K e^{\beta_k x_i}}}
#' @return An object of class `params_mlogit`, which is a list containing `coefs`
#' and `n_samples`, where `n_samples` is equal to the number of rows
#' in each element of `coefs`.
#' @examples 
#' params <- params_mlogit(
#'   coefs = list(
#'     healthy_to_sick = data.frame(
#'       intercept = c(-0.33, -.2, -.15),
#'       treat = c(log(.75), log(.8), log(.9))
#'     ),
#'    
#'     healthy_to_death = data.frame(
#'       intercept = c(-1, -1.2, -.5),
#'       treat = c(log(.6), log(.65), log(.55))
#'     )
#'   )
#' )
#' summary(params)
#' params
#' @aliases print.params_mlogit
#' @export
params_mlogit <- function(coefs){
  if (inherits(coefs, "list")) coefs <- coeflist_to_array(coefs)
  if (is.null(colnames(coefs))) colnames(coefs) <- paste0("x", 1:ncol(coefs))
  n_samples <- get_n_samples(coefs)
  return(check(new_params_mlogit(coefs, n_samples)))
}

new_params_mlogit <- function(coefs, n_samples){
  stopifnot(is.numeric(n_samples))
  res <- list(coefs = coefs, n_samples = n_samples)
  class(res) <- "params_mlogit"
  return(res)
}

#' @rdname check
check.params_mlogit <- function(object){
  if (!(is.numeric(object$coefs) && is_3d_array(object$coefs))) {
    stop("'coefs' must be a numeric 3D array.")
  } 
  return(object)
}

# summary.params_mlogit() --------------------------------------------------------
#' @export
#' @rdname summary.params
summary.params_mlogit <- function(object, prob = 0.95, ...) {
  
 rbindlist(apply(object$coef, 3, coef_summary, prob = prob),
          idcol = "transition")
}

# print.params_mlogit() --------------------------------------------------------
#' @export
print.params_mlogit <- function(x, ...) {
  
  cat("A \"params_mlogit\" object:\n\n")
  cat("Summary of coefficient estimates:\n")
  print(summary(x))
  cat("\n")
  cat(paste0("Number of parameter samples: ", x$n_samples))
  cat("\n")
  cat(paste0("Number of health states: ",  dim(x$coefs)[3] + 1))

  invisible(x)
}

# create_params.multinom() -----------------------------------------------------
create_coef_multinom <- function(object, n = 1000, uncertainty){
  # Extract/simulate coefficients
  coefs <- c(t(stats::coef(object)))
  if (uncertainty == "normal"){
    coefs_sim <- MASS::mvrnorm(n = n, 
                               mu = coefs,
                               Sigma = stats::vcov(object))
    if (n == 1) coefs_sim <- matrix(coefs_sim, nrow = 1)
  } else{
    coefs_sim <- matrix(coefs, nrow = 1)
  }
  coefs_est <- stats::coef(object) # The point estimate
  if (is_1d_vector(coefs_est)){
    coefs_est <- matrix(coefs_est, nrow = 1)
    colnames(coefs_est) <- names(stats::coef(object))
  } 
  
  # Store coefficients in array
  p <- ncol(coefs_est)
  K <- nrow(coefs_est)
  coefs_sim_array <- array(NA, dim = c(nrow(coefs_sim), p, K))
  start_col <- 1
  end_col <- p
  for (j in 1:K){
    coefs_sim_array[, , j] <- coefs_sim[, start_col:end_col, drop = FALSE]
    start_col <- end_col + 1
    end_col <- end_col + p
  } 
  
  # Array dimension names
  class_names <- NULL
  if(length(object$lev)) class_names < object$lev[-1L]
  if(length(object$lab)) class_names <- object$lab[-1L]
  dimnames(coefs_sim_array) <- list(NULL,
                                    colnames(coefs_est),
                                    class_names)
  
  # Return
  return(coefs_sim_array)
}

#' @export
#' @rdname create_params
create_params.multinom <- function(object, n = 1000, uncertainty = c("normal", "none"),
                                   ...){
  uncertainty <- deprecate_point_estimate(list(...)$point_estimate, uncertainty,
                                          missing(uncertainty))
  uncertainty <- match.arg(uncertainty)
  coefs <- create_coef_multinom(object, n, uncertainty)
  if (uncertainty == "none"){
    n_samples <- 1
  } else{
    n_samples <- n
  }
  return(new_params_mlogit(coefs = coefs,
                           n_samples = n_samples))
}