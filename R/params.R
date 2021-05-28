# Generic documentation for parameter object -----------------------------------
#' Parameter object
#' 
#' Objects prefixed by "params_" are lists containing the parameters of a statistical model
#' used for simulation modeling. The parameters are used to simulate outcomes
#' as a function of covariates contained in input matrices ([input_mats]).
#' 
#' @name params
#' @seealso [`tparams`]
NULL

# Helper functions -------------------------------------------------------------
new_params_list <- function(..., inner_class, new_class){
  return(object_list(..., inner_class = inner_class,
                     new_class = new_class))
}

check_params_list <- function(x){
  inner_class <- class(x[[1]])
  n_samples <- sapply(x, function(y) y$n_samples)
  if(!all(n_samples == n_samples[1])){
    msg <- paste0("The number of samples in each '", inner_class, "' object must be the same.")
    stop(msg, call. = FALSE)
  }
  return(x)
}

create_params_list <- function(object, n, uncertainty, inner_class, new_class,
                               ...){
  n_objects <- length(object)
  params_list <- vector(mode = "list", length = n_objects)
  names(params_list) <- names(object)
  for (i in 1:n_objects){
    params_list[[i]] <- create_params(object[[i]], n, uncertainty, ...)
  }
  return(new_params_list(params_list, inner_class = inner_class,
                         new_class = new_class))
}

coef_summary <- function(x, prob) {
  alpha <- ci_alpha(prob)
  data.table(
    term = colnames(x),
    estimate = apply(x, 2, mean),
    lower = apply(x, 2, stats::quantile, prob = alpha$lower),
    upper = apply(x, 2, stats::quantile, prob = alpha$upper)
  )
}

summary_params_list <- function(object, prob = 0.95, idcol = "model", ...) {
  res <- lapply(object, summary)
  rbindlist(res, idcol = idcol)
}

# create_params generic method -------------------------------------------------
#' Create a parameter object from a fitted model
#' 
#' `create_params` is a generic function for creating an object containing 
#' parameters from a fitted statistical model. If `uncertainty != "none"`,
#' then random samples from suitable probability distributions are returned.
#' @param object A statistical model to randomly sample parameters from.  
#' @param n Number of random observations to draw. Not used if `uncertainty = "none"`.
#' @param uncertainty Method determining how parameter uncertainty should be handled. 
#' If `"normal"`, then parameters are randomly drawn from their multivariate normal
#' distribution. If `"bootstrap`, then parameters are bootstrapped using [`bootstrap`].
#' If `"none`, then only point estimates are returned.
#' @param ... Further arguments passed to or from other methods. Only used when
#'  `object` is of class `partsurvfit`, in which case the arguments are passed 
#'  to [`bootstrap.partsurvfit()`].
#' @return An object prefixed by `params_`. Mapping between `create_params` 
#' and the classes of the returned objects are: 
#' \itemize{
#' \item{`create_params.lm` ->}{ [`params_lm`]}
#' \item{`create_params.multinom` ->}{ [`params_mlogit`]}
#' \item{`create_params.multinom_list` ->}{ [`params_mlogit_list`]}
#' \item{`create_params.flexsurvreg` ->}{ [`params_surv`]}
#' \item{`create_params.flexsurvreg_list` ->}{ [`params_surv_list`]}
#' \item{`create_params.partsurvfit` ->}{ [`params_surv_list`]}
#' }
#' @name create_params
#' @examples 
#' # create_params.lm
#' fit <- stats::lm(costs ~ female, data = psm4_exdata$costs$medical)
#' n <- 5
#' params_lm <- create_params(fit, n = n)
#' head(params_lm$coefs)
#' head(params_lm$sigma)
#' 
#' # create_params.flexsurvreg
#' library("flexsurv")
#' fit <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                     data = ovarian, dist = "weibull")
#' n <- 5
#' params_surv_wei <- create_params(fit, n = n)
#' print(params_surv_wei$dist)
#' head(params_surv_wei$coefs)
#' @export
#' @rdname create_params
create_params <- function (object, ...) {
  UseMethod("create_params", object)
}

# Linear model -----------------------------------------------------------------
#' Parameters of a linear model
#' 
#' Create a list containing the parameters of a fitted linear regression model.
#' @param coefs  Matrix of samples of the coefficients under sampling uncertainty.
#' @param sigma A vector of samples of the standard error of the regression model. 
#' Must only be specified if the model is used to randomly simulate values 
#' (rather than to predict means).
#' 
#' @details Fitted linear models are used to predict values, \eqn{y},
#'  as a function of covariates, \eqn{x},
#' \deqn{y = x^T\beta + \epsilon.}
#' Predicted means are given by \eqn{x^T\hat{\beta}} where \eqn{\hat{\beta}}
#' is the vector of estimated regression coefficients. Random samples are obtained by 
#' sampling the error term from a normal distribution, 
#' \eqn{\epsilon \sim N(0, \hat{\sigma}^2)}{\epsilon ~ N(0, \hat{\sigma}^2)}.
#' @return An object of class `params_lm`, which is a list containing `coefs`,
#' `sigma`, and `n_samples`. `n_samples` is equal to the number of rows
#' in `coefs`.
#' @examples 
#' library("MASS")
#' n <- 2
#' params <- params_lm(coefs = MASS::mvrnorm(n, mu = c(.5,.6),
#'                                             Sigma = matrix(c(.05, .01, .01, .05), nrow = 2)),
#'                       sigma <- rgamma(n, shape = .5, rate = 4))
#' print(params)
#'
#' @export
params_lm <- function(coefs, sigma = NULL){
  stopifnot(is.matrix(coefs) | is.vector(coefs))
  if (is.vector(coefs)) coefs <- matrix(coefs, ncol = 1)
  if(is.null(sigma)) sigma <- rep(1, nrow(coefs))
  n_samples <- nrow(coefs)
  check(new_params_lm(coefs, sigma, n_samples))
}

new_params_lm <- function(coefs, sigma, n_samples){
  stopifnot(is.numeric(sigma))
  stopifnot(is.numeric(n_samples))
  l <- list(coefs = coefs, sigma = sigma, n_samples = n_samples)
  class(l) <- "params_lm"
  return(l)
}

#' @rdname check
check.params_lm <- function(object){
  if(object$n_samples != length(object$sigma)){
    stop("Number of samples in 'sigma' is not equal to the number of samples in 'coefs'.",
         call. = FALSE)
  }
  return(object)
}

#' Parameters of a list of linear models
#' 
#' Create a list containing the parameters of a list of fitted linear regression models.
#' @param ... Objects of class [`params_lm`], which can be named.
#' 
#' @return An object of class `params_lm_list`, which is a list containing [`params_lm`]
#' objects. 
#' @export
#' @keywords internal
#' @examples 
#' n <- 2
#' coefs1 <- MASS::mvrnorm(n, mu = c(.5,.6),
#'                         Sigma = matrix(c(.05, .01, .01, .05), nrow = 2))
#' sigma1 <- rgamma(n, shape = .5, rate = 4)
#' params1 <- params_lm(coefs = coefs1, sigma = sigma1)

#' coefs2 <- MASS::mvrnorm(n, mu = c(.2,.9),
#'                         Sigma = matrix(c(.08, .02, .02, .08), nrow = 2))
#' sigma2 <- rgamma(n, shape = .9, rate = 4)
#' params2 <- params_lm(coefs = coefs2, sigma = sigma2)

#' params_list <- params_lm_list(params1, params2)
#' print(params_list)
params_lm_list <- function(...){
  return(check(new_params_list(..., inner_class = "params_lm",
                               new_class = "params_lm_list")))
}

#' @rdname check
check.params_lm_list <- function(object){
  check_params_list(object)
}

#' @export
#' @rdname create_params
create_params.lm <- function(object, n = 1000, uncertainty = c("normal", "none"),
                             ...){
  uncertainty <- deprecate_point_estimate(list(...)$point_estimate, uncertainty,
                                          missing(uncertainty))
  uncertainty <- match.arg(uncertainty)
  
  if (uncertainty == "normal"){
    coefs_sim <- MASS::mvrnorm(n, stats::coef(object), stats::vcov(object))
    if(is.vector(coefs_sim)) coefs_sim <- matrix(coefs_sim, nrow = 1)
    return(new_params_lm(coefs = coefs_sim,
                         sigma = rep(summary(object)$sigma, n),
                         n_samples = n))
  } else{
    coefs <- matrix(stats::coef(object), nrow = 1)
    colnames(coefs) <- names(stats::coef(object))
    return(new_params_lm(coefs = coefs,
                         sigma = summary(object)$sigma,
                         n_samples = 1))
  }
}

#' @export
# #' @rdname create_params
#' @keywords internal
create_params.lm_list <- function(object, n = 1000, uncertainty = c("normal", "none"),
                                   ...){
  return(create_params_list(object, n = n, uncertainty = uncertainty, 
                            inner_class = "params_lm", new_class = "params_lm_list"))
}

