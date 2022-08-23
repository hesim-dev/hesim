# params_lm() ------------------------------------------------------------------
#' Parameters of a linear model
#'
#' Create a list containing the parameters of a fitted linear regression model.
#' @param coefs  Samples of the coefficients under sampling uncertainty.
#' Must be a matrix or any object coercible to a matrix such as `data.frame`
#' or `data.table`.
#' @param sigma A vector of samples of the standard error of the regression model.
#' Default value is 1 for all samples. Only used if the model is
#' used to randomly simulate values (rather than to predict means).
#'
#' @details Fitted linear models are used to predict values, \eqn{y},
#'  as a function of covariates, \eqn{x},
#' \deqn{y = x^T\beta + \epsilon.}
#' Predicted means are given by \eqn{x^T\hat{\beta}} where \eqn{\hat{\beta}}
#' is the vector of estimated regression coefficients. Random samples are obtained by
#' sampling the error term from a normal distribution,
#' \eqn{\epsilon \sim N(0, \hat{\sigma}^2)}{\epsilon ~ N(0, \hat{\sigma}^2)}.
#' @return An object of class `params_lm`, which is a list containing `coefs`,
#' `sigma`, and `n_samples`. `n_samples` is equal to the
#'  number of rows in `coefs`. The `coefs` element is always converted into a
#'  matrix.
#' @examples
#' library("MASS")
#' n <- 2
#' params <- params_lm(
#'   coefs = mvrnorm(n,
#'     mu = c(.5, .6),
#'     Sigma = matrix(c(.05, .01, .01, .05), nrow = 2)
#'   ),
#'   sigma <- rgamma(n, shape = .5, rate = 4)
#' )
#' summary(params)
#' params
#'
#' @seealso This parameter object is useful for modeling [health state values][StateVals]
#' when values can vary across patients and/or health states as a function of
#' covariates. In many cases it will, however, be simpler, and more flexible to
#' use a [`stateval_tbl`]. For an example use case see the documentation for
#' [create_StateVals.lm()].
#' @export
params_lm <- function(coefs, sigma = 1) {
  coefs <- as.matrix(coefs)
  if (length(sigma) == 1) sigma <- rep(sigma, nrow(coefs))
  n_samples <- nrow(coefs)
  check(new_params_lm(coefs, sigma, n_samples))
}

new_params_lm <- function(coefs, sigma, n_samples) {
  stopifnot(is.numeric(sigma))
  stopifnot(is.numeric(n_samples))
  if (is.null(colnames(coefs))) colnames(coefs) <- default_colnames(coefs)
  l <- list(coefs = coefs, sigma = sigma, n_samples = n_samples)
  class(l) <- "params_lm"
  return(l)
}

#' @rdname check
check.params_lm <- function(object) {
  if (object$n_samples != length(object$sigma)) {
    stop("Number of samples in 'sigma' is not equal to the number of samples in 'coefs'.",
      call. = FALSE
    )
  }
  return(object)
}

# summary.params_lm() ----------------------------------------------------------
#' @rdname summary.params
#' @export
summary.params_lm <- function(object, probs = c(.025, .975), ...) {
  sigma_mat <- matrix(object$sigma, ncol = 1)
  colnames(sigma_mat) <- "sigma"

  rbindlist(
    lapply(list(mean = object$coef, sd = sigma_mat),
      coef_summary,
      probs = probs
    ),
    idcol = "parameter"
  )
}

# print.params_lm() ------------------------------------------------------------
#' @export
print.params_lm <- function(x, ...) {
  parameter <- NULL
  x_summary <- summary(x)

  cat("A \"params_lm\" object\n\n")
  cat("Summary of coefficients:\n")
  print(x_summary[parameter == "mean"])
  cat("\n")
  cat("Summary of sigma:\n")
  print(x_summary[parameter == "sd"])
  invisible(x)
}

# create_params.lm() -----------------------------------------------------------
#' @export
#' @rdname create_params
create_params.lm <- function(object, n = 1000, uncertainty = c("normal", "none"),
                             ...) {
  uncertainty <- deprecate_point_estimate(
    list(...)$point_estimate, uncertainty,
    missing(uncertainty)
  )
  uncertainty <- match.arg(uncertainty)

  if (uncertainty == "normal") {
    coefs_sim <- MASS::mvrnorm(n, stats::coef(object), stats::vcov(object))
    if (n == 1) {
      coefs_sim <- matrix(coefs_sim, nrow = 1)
      colnames(coefs_sim) <- names(stats::coef(object))
    }
    return(new_params_lm(
      coefs = coefs_sim,
      sigma = rep(summary(object)$sigma, n),
      n_samples = n
    ))
  } else {
    coefs <- matrix(stats::coef(object), nrow = 1)
    colnames(coefs) <- names(stats::coef(object))
    return(new_params_lm(
      coefs = coefs,
      sigma = summary(object)$sigma,
      n_samples = 1
    ))
  }
}
