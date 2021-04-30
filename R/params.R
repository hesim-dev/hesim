# Generic parameter object documentation ---------------------------------------
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

new_params_joined <- function(..., times, inner_class){
  objects <- create_object_list(...)
  res <- list(models = objects, times = times)
  class(res) <- paste0("params_joined", sub("params", "", inner_class))
  return(res)
}

check_params_joined <- function(x, inner_class, model_list){
  check_joined_object(x, inner_class = inner_class, model_list = model_list)
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

create_params_joined <- function(object, n, uncertainty, inner_class){
  n_models <- length(object$models)
  models <- vector(mode = "list", length = n_models)
  names(models) <- names(object$models)
  for (i in 1:n_models){
    models[[i]] <- create_params(object$models[[i]], n = n, uncertainty = uncertainty)
  }
  return(new_params_joined(models, times = object$times, inner_class = inner_class))
}

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


# Multinomial logit model ------------------------------------------------------
#' Parameters of a multinomial logit model
#' 
#' Store the parameters of a fitted multinomial logistic 
#' regression model. The model is used to predict probabilities of \eqn{K} 
#' classes.
#' @param coefs  A 3D array of stacked matrices. The number of matrices (i.e.,
#' the number of slices in the cube) should be equal to \eqn{K-1}. Each 
#' matrix is contains samples of the regression coefficients under sampling uncertainty
#' corresponding to a particular class. Rows index parameter samples and 
#'  columns index coefficients.
#' 
#' @details Multinomial logit models are used to predict the probability of 
#' membership for subject \eqn{i} in each of \eqn{K} classes as a function of covariates:
#' \deqn{Pr(y_i = c) = \frac{e^{\beta_c x_i}}{\sum_{k=1}^K e^{\beta_k x_i}}}
#' @return An object of class `params_mlogit`, which is a list containing `coefs`
#' and `n_samples`, where `n_samples` is equal to the number of rows
#' in each element of `coefs`.
#' @examples 
#' params <- params_mlogit(coefs = array(
#'   c(matrix(c(intercept = 0, treatment = log(.75)), nrow = 1),
#'     matrix(c(intercept = 0, treatment = log(.8)), nrow = 1)),
#'   dim = c(1, 2, 2)
#' )) 
#' @export
params_mlogit <- function(coefs){
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
  if (!is.numeric(object$coefs)) stop("'coefs' must be a numeric array.")
  return(object)
}

#' Parameters of a list of multinomial logit models
#' 
#' Create a list containing the parameters of multiple fitted multinomial logit models.
#' @param ... Objects of class [`params_mlogit`], which can be named.
#' 
#' @return An object of class `params_mlogit_list`, which is a list containing 
#' [`params_mlogit`] objects.
#' @export
params_mlogit_list <- function(...){
  return(check(new_params_list(..., inner_class = "params_mlogit", 
                               new_class = "params_mlogit_list")))
}

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

#' @export
#' @rdname create_params
create_params.multinom_list <- function(object, n = 1000, uncertainty = c("normal", "none"), ...){
  return(create_params_list(object, n = n, uncertainty = uncertainty, 
                            inner_class = "params_mlogit", new_class = "params_mlogit_list",
                            ...))
}

# List of survival models ------------------------------------------------------
#' Parameters of a list of survival models
#' 
#' Create a list containing the parameters of multiple fitted parametric survival models.
#' @param ... Objects of class [`params_surv`], which can be named.
#' 
#' @return An object of class `params_surv_list`, which is a list containing [`params_surv`]
#' objects.
#' @examples 
#' library("flexsurv")
#' fit_wei <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist = "weibull")
#' params_wei <- create_params(fit_wei, n = 2)
#' 
#' fit_exp <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist = "exp")
#' params_exp <- create_params(fit_exp, n = 2)
#' 
#' params_list <- params_surv_list(wei = params_wei, exp = params_exp)
#' print(params_list)
#' @export
params_surv_list <- function(...){
  return(check(new_params_list(..., inner_class = "params_surv", 
                               new_class = "params_surv_list")))
}

#' @rdname check
check.params_surv_list <- function(object){
  check_params_list(object)
}

#' @export
#' @rdname create_params
create_params.flexsurvreg_list <- function(object, n = 1000, uncertainty = c("normal", "none"),
                                           ...){
  return(create_params_list(object, n = n, uncertainty = uncertainty, 
                            inner_class = "params_surv", new_class = "params_surv_list"))
}

#' @export
#' @rdname create_params
create_params.partsurvfit <- function(object, n = 1000, 
                                      uncertainty = c("normal", "bootstrap", "none"), 
                                      ...){
  uncertainty <- match.arg(uncertainty)
  if(uncertainty == "bootstrap"){
    res <- bootstrap(object, B = n, ...)
  } else{
    res <- create_params(object$models, n = n, uncertainty = uncertainty)
  }
  return(res)
}

# Joined survival models -------------------------------------------------------
#' Parameters of joined survival models
#' 
#' Create a list containing the parameters of survival models joined at specified time points. See
#' \code{\link{joined}} for more details.
#' @param ... Objects of class [`params_surv`], which can be named.
#' @param times A numeric vector of times at which to join models.
#' 
#' @return An object of class `params_joined_surv`, which is a list containing two elements:
#' \describe{
#' \item{models}{A list of [`params_surv`] objects from each statistical model
#' to be joined.}
#' \item{times}{Equivalent to the argument `times`.}
#' }
#' @examples 
#' library("flexsurv")
#' fit_exp <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                                  data = ovarian, dist = "exp")
#' fit_wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                                  data = ovarian, dist = "weibull")
#' params_surv_exp <- create_params(fit_exp, n = 2)
#' params_surv_wei <- create_params(fit_wei, n = 2)
#' params_joined_surv <- params_joined_surv(exp = params_surv_exp,
#'                                          wei = params_surv_wei,
#'                                          times = 3)
#' print(params_joined_surv)
#' @export
params_joined_surv <- function(..., times){
  return(check(new_params_joined(..., times = times, inner_class = "params_surv"),
               inner_class = "params_surv"))
}

#' @rdname check
check.params_joined_surv <- function(object, inner_class){
  check_params_joined(object, inner_class = inner_class, model_list = FALSE)
}

# Joined list of survival models -----------------------------------------------
#' Parameters of joined lists of survival models
#' 
#' Create a list containing the parameters of multiple sets of survival models, each joined
#' at specified time points. See [`joined`] for more details.
#' @param ... Objects of class [`params_surv_list`], which can be named.
#' @param times A list of sorted numeric vectors, with the length of each list element equal
#' to the number of sets of models.
#' 
#' @return An object of class `params_joined_surv_list`, which is a list containing two elements:
#' \describe{
#' \item{models}{A list of [`params_surv_list`], each containing [`params_surv`] 
#' objects to be joined.}
#' \item{times}{Equivalent to the argument `times`.}
#' }
#' @examples 
#' library("flexsurv")
#' fit_exp <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                                  data = ovarian, dist = "exp")
#' fit_wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                                  data = ovarian, dist = "weibull")
#' fit_lnorm <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                                   data = ovarian, dist = "lognormal")
#'
#' params_exp <- create_params(fit_exp, n = 2)
#' params_wei <- create_params(fit_wei, n = 2)
#' params_lnorm <- create_params(fit_lnorm, n = 2)
#'
#' params_list1 <- params_surv_list(params_exp, params_wei)
#' params_list2 <- params_surv_list(params_exp, params_lnorm)
#' params_joined <- params_joined_surv_list(model1 = params_list1,
#'                                          model2 = params_list2,
#'                                          times = list(3, 5))
#' print(params_joined)
#' @export
params_joined_surv_list <- function(..., times){
  return(check(new_params_joined(..., times = times, inner_class = "params_surv_list"),
               inner_class = "params_surv_list"))
}

#' @rdname check
check.params_joined_surv_list <- function(object, inner_class){
  check_params_joined(object, inner_class = inner_class, model_list = TRUE)
}

#' @export
# #' @rdname create_params
#' @keywords internal
create_params.joined_flexsurvreg <- function(object, n = 1000, uncertainty = c("normal", "none"),
                                             ...){
  return(create_params_joined(object, n = n, uncertainty = uncertainty,
                              inner_class =  "params_surv"))
}

#' @export
# #' @rdname create_params
#' @keywords internal
create_params.joined_flexsurvreg_list <- function(object, n = 1000, uncertainty = c("normal", "none"),
                                                  ...){
  return(create_params_joined(object, n = n, uncertainty = uncertainty, 
                              inner_class = "params_surv_list"))
}

