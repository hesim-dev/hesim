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

#' Parameters of a mean model
#' 
#' Create a list containing the parameters of a mean model.
#' @param mu  Matrix of samples from the posterior distribution of the 
#' mean. Columns denote random samples and rows denote means for different observations.
#' @param sigma A vector of samples of the standard deviation.
#' 
#' @return An object of class "params_mean", which is a list containing \code{mu},
#' \code{sigma}, and \code{n_samples}.
#'  \code{n_samples} is equal to the number of columns in \code{mu}.
#' @examples 
#' params <- params_mean(mu = matrix(seq(1, 4), nrow = 2), 
#'                       sigma = c(0, 0))
#' print(params)
#'
#' @export
params_mean <- function(mu, sigma = NULL){
  stopifnot(is.matrix(mu))
  n_samples <- ncol(mu)
  if(is.null(sigma)) sigma <- rep(1, n_samples)
  check(new_params_mean(mu, sigma, n_samples))
}

new_params_mean <- function(mu, sigma, n_samples, strategy_id, patient_id, state_id){
  stopifnot(is.numeric(sigma))
  stopifnot(is.numeric(n_samples))
  l <- list(mu = mu, sigma = sigma, n_samples = n_samples)
  class(l) <- "params_mean"
  return(l)
}

#' @rdname check
check.params_mean <- function(object){
  if(object$n_samples != length(object$sigma)){
    stop("Number of samples in 'sigma' is not equal to the number of samples in 'mu'.",
         call. = FALSE)
  }
  return(object)
}

#' Parameters of a linear model
#' 
#' Create a list containing the parameters of a fitted linear regression model.
#' @param coefs  Matrix of samples from the posterior distribution of the 
#' regression coefficients.
#' @param sigma A vector of samples of the standard error of the regression model.
#' 
#' @return An object of class "params_lm", which is a list containing \code{coefs},
#' \code{sigma}, and \code{n_samples}. \code{n_samples} is equal to the number of rows
#' in \code{coefs}.
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
#' @param ... Objects of class \code{\link{params_lm}}, which can be named.
#' 
#' @return An object of class "params_lm_list", which is a list containing \code{\link{params_lm}}
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

#' Parameters of a survival model
#' 
#' Create a list containing the parameters of a single fitted parametric survival model.
#' @param coefs A list of length equal to the number of parameters in the 
#' survival distribution. Each element of the list is a matrix of samples from 
#' the posterior distribution of the regression coefficients used to predict
#' a given parameter.
#' @param dist Character vector denoting the parametric distribution. See "Details".
#' @param aux Auxillary arguments used with splines or fractional polynomials. See "Details". 
#' 
#' @return An object of class "params_surv", which is a list containing \code{coefs},
#' \code{dist}, and \code{n_samples}. \code{n_samples} is equal to the number of rows
#' in each element of \code{coefs}, which must be the same. The list may also contain \code{aux} if
#' a spline or fractional polynomial model is fit. 
#' 
#' @details 
#' The types of distributions that can be specified are: 
#' \itemize{
#' \item{\code{exponential} or \code{exp}}{ Exponential distribution. \code{coef}
#' must contain the \code{rate} parameter on the log scale and the same parameterization as in 
#' \code{\link[stats]{Exponential}}}.
#' \item{\code{weibull} or \code{weibull.quiet}}{ Weibull distribution. The first 
#' element of \code{coef} is the \code{shape} parameter (on the log scale) and the second
#' element is the \code{scale} parameter (also on the log scale). The parameterization is
#' that same as in \code{\link[stats]{Weibull}}.}
#' \item{\code{gamma}}{ Gamma distribution. The first 
#' element of \code{coef} is the \code{shape} parameter (on the log scale) and the second
#' element is the \code{rate} parameter (also on the log scale). The parameterization is
#' that same as in \code{\link[stats]{GammaDist}}.}
#' \item{\code{lnorm}}{ Lognormal distribution. The first 
#' element of \code{coef} is the \code{meanlog} parameter (i.e., the mean on the log scale) and the second
#' element is the \code{sdlog} parameter (i.e., the standard deviation on the log scale). The parameterization is
#' that same as in \code{\link[stats]{Lognormal}}.}
#' \item{\code{gompertz}}{ Gompertz distribution. The first 
#' element of \code{coef} is the \code{shape} parameter and the second
#' element is the \code{rate} parameter (on the log scale). The parameterization is
#' that same as in \code{\link[flexsurv]{Gompertz}}.}
#' \item{\code{llogis}}{ Log-logistic distribution. The first 
#' element of \code{coef} is the \code{shape} parameter (on the log scale) and the second
#' element is the \code{scale} parameter (also on the log scale). The parameterization is
#' that same as in \code{\link[flexsurv]{Llogis}}.}
#' \item{\code{gengamma}}{ Generalized gamma distribution. The first 
#' element of \code{coef} is the location parameter \code{mu}, the second
#' element is the scale parameter \code{sigma} (on the log scale), and the
#' third element is the shape parameter \code{Q}. The parameterization is
#' that same as in \code{\link[flexsurv]{GenGamma}}.}
#' \item{\code{survspline}}{ Survival splines. Each element of \code{coef} is a parameter of the
#' spline model (i.e. \code{gamma_0}, \code{gamma_1}, \eqn{\ldots}) with length equal
#' to the number of knots (including the boundary knots). See below for details on the
#' auxillary arguments. The parameterization is that same as in \code{\link[flexsurv]{Survspline}}.}
#' \item{\code{fracpoly}}{ Fractional polynomials. Each element of \code{coef} is a parameter of the
#' fractional polynomial model (i.e. \code{gamma_0}, \code{gamma_1}, \eqn{\ldots}) with length equal
#' to the number of powers minus 1. See below for details on the auxillary arguments 
#' (i.e., \code{powers}).}
#' }
#' 
#' Auxillary arguments for spline models should be specified as a list containing the elements:
#' \describe{
#' \item{\code{knots}}{A numeric vector of knots.}
#' \item{\code{scale}}{The survival outcome to be modeled
#' as a spline function. Options are "log_cumhazard" for the log cummulative hazard; 
#' "log_hazard" for the log hazard rate; "log_cumodds" for the log cummulative odds;
#' and "inv_normal" for the inverse normal distribution function.}
#' \item{\code{timescale}}{If "log" (the default), then survival is modeled as a spline function
#' of log time; if "identity", then it is modeled as a spline function of time.}
#' }
#' 
#' Auxillary arguments for fractional polynomial models should be specified as a list containing the elements:
#' \describe{
#' \item{\code{powers}}{ A vector of the powers of the fractional polynomial with each element
#'  chosen from the following set: -2. -1, -0.5, 0, 0.5, 1, 2, 3.}
#' }
#' 
#' Furthermore, when either splines for fractional polynomials are used, the following additional auxillary arguments
#' can be specified:
#' \describe{
#' \item{\code{cumhaz_method}}{Numerical method used to compute cumulative hazard 
#' (i.e., to integrate the hazard function). Options are and "quad" for adaptive
#' quadrature and "riemann" for Riemann sum.}
#' \item{\code{step}}{Step size for computation of cumulative hazard with 
#' numerical integration. Only required when integrating with
#'  using "riemann". No step size is required for "quad".}
#' }
#' 
#' @examples 
#' library("flexsurv")
#' fit <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist = "weibull")
#' params <- params_surv(coefs = list(shape = fit$res.t["shape", "est", drop = FALSE],
#'                                    scale = fit$res.t["scale", "est", drop = FALSE]),
#'                      dist = fit$dlist$name)
#' print(params)
#' @export
params_surv <- function(coefs, dist, aux = NULL){
  stopifnot(is.list(coefs))
  if (!is.matrix(coefs[[1]])){
    stop("'coefs' must be a list of matrices.", call. = FALSE)
  }
  n_samples <- nrow(coefs[[1]])
  check(new_params_surv(coefs, dist, n_samples, aux))
}

new_params_surv <- function(coefs, dist, n_samples, aux = NULL){
  stopifnot(is.vector(dist) & is.character(dist))
  stopifnot(is.numeric(n_samples))
  res <- list(coefs = coefs, dist = dist)
  if (!is.null(aux)) {
    res[["aux"]] <- aux
    if (is.null(aux$cumhaz_method)){
      res[["aux"]]$cumhaz_method <- "quad"
    }
    if (is.null(aux$step)){
      if(res[["aux"]]$cumhaz_method %in% c("riemann")){
        msg <- paste0("If the Riemann sum is used to compute the cumulative ",
                     "hazard, then the step size must be specified.")
        stop(msg)
      } 
    }    
  }
  res[["n_samples"]] <- n_samples
  class(res) <- "params_surv"
  return(res)
}

#' @rdname check
check.params_surv <- function(object){
  if (list_depth(object$coefs) !=1 | length(object$dist) !=1){
    stop("'coefs' must only contain one survival model.", call. = FALSE)
  }
  matrix_bool <- unlist(lapply(object$coefs, is.matrix))
  if(sum(!matrix_bool) > 0){
    stop("'coefs' must be a list of matrices.",
         call. = FALSE)
  }
  coefs_nrows <- unlist(lapply(object$coefs, nrow))
  if(!all(object$n_samples == coefs_nrows)){
    stop("Number of rows in all 'coefs' matrices must be equal.",
         call. = FALSE)
  } 
  return(object)
}

#' Parameters of joined survival models
#' 
#' Create a list containing the parameters of survival models joined at specified time points. See
#' \code{\link{joined}} for more details.
#' @param ... Objects of class \code{\link{params_surv}}, which can be named.
#' @param times A numeric vector of times at which to join models.
#' 
#' @return An object of class "params_joined_surv", which is a list containing two elements:
#' \describe{
#' \item{models}{A list of \code{\link{params_surv}} objects from each statistical model
#' to be joined.}
#' \item{times}{Equivalent to the argument \code{times}.}
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

#' Parameters of a list of survival models
#' 
#' Create a list containing the parameters of multiple fitted parametric survival models.
#' @param ... Objects of class \code{\link{params_surv}}, which can be named.
#' 
#' @return An object of class "params_surv_list", which is a list containing \code{params_surv}
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

#' Parameters of joined lists of survival models
#' 
#' Create a list containing the parameters of multiple sets of survival models, each joined
#' at specified time points. See \code{\link{joined}} for more details.
#' @param ... Objects of class \code{\link{params_surv_list}}, which can be named.
#' @param times A list of sorted numeric vectors, with the length of each list element equal
#' to the number of sets of models.
#' 
#' @return An object of class "params_joined_surv_list", which is a list containing two elements:
#' \describe{
#' \item{models}{A list of \code{\link{params_surv_list}}, each containing code{\link{params_surv}} 
#' objects to be joined.}
#' \item{times}{Equivalent to the argument \code{times}.}
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

create_params_list <- function(object, n, point_estimate, inner_class, new_class){
  n_objects <- length(object)
  params_list <- vector(mode = "list", length = n_objects)
  names(params_list) <- names(object)
  for (i in 1:n_objects){
    params_list[[i]] <- create_params(object[[i]], n, point_estimate)
  }
  return(new_params_list(params_list, inner_class = inner_class,
                         new_class = new_class))
}

create_params_joined <- function(object, n, point_estimate, inner_class){
  n_models <- length(object$models)
  models <- vector(mode = "list", length = n_models)
  names(models) <- names(object$models)
  for (i in 1:n_models){
    models[[i]] <- create_params(object$models[[i]], n, point_estimate)
  }
  return(new_params_joined(models, times = object$times, inner_class = inner_class))
}

#' Create a parameter object from a fitted model
#' 
#' \code{create_params} is a generic function for creating an object containing 
#' parameters from a fitted statistical model. If \code{point_estimate = FALSE},
#' then random samples from the posterior distribution are returned.
#' @param object A statistical model to randomly sample parameters from.  
#' @param n Number of random observations to draw. 
#' @param point_estimate If TRUE, then the point estimates are returned and
#' and no samples are drawn.
#' @param bootstrap If \code{bootstrap} is FALSE or not specified, then \code{n} parameter sets are 
#' drawn by sampling from a multivariate normal distribution. If \code{bootstrap} is TRUE, then 
#' parameters are bootstrapped using \code{\link{bootstrap}}. 
#' @param max_errors Equivalent to the \code{max_errors} argument in \code{\link{bootstrap}}. 
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' @return An object prefixed by \code{params_}. Mapping between \code{create_params} 
#' and the classes of the returned objects are: 
#' \itemize{
#' \item{\code{create_params.lm} ->}{ \code{params_lm}}
#' \item{\code{create_params.flexsurvreg} ->}{ \code{params_surv}}
#' \item{\code{create_params.flexsurvreg_list} ->}{ \code{params_surv_list}}
#' \item{\code{create_params.partsurvfit} ->}{ \code{params_surv_list}}
#' }
#' @keywords internal
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

#' @export
#' @rdname create_params
create_params.lm <- function(object, n = 1000, point_estimate = FALSE, ...){
  if (point_estimate == FALSE){
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
create_params.lm_list <- function(object, n = 1000, point_estimate = FALSE, ...){
  return(create_params_list(object, n, point_estimate, 
                          inner_class = "params_lm", new_class = "params_lm_list"))
}

flexsurvreg_inds <- function(object){
  n_pars <- length(object$dlist$pars)
  inds <- vector(mode = "list", length = n_pars)
  for (j in 1:n_pars){
    parname_j <-  object$dlist$pars[j]
    covind_j <- object$mx[[parname_j]]
    if (length(covind_j) > 0){
      inds[[j]] <- c(j, n_pars + covind_j)
    } else{
      inds[[j]] <- j
    }
  }
  return(inds)
}

flexsurvreg_aux <- function(object){
  if(object$dlist$name == "survspline"){
    aux <- object$aux
    if(aux$scale == "hazard") aux$scale <- "log_cumhazard"
    if(aux$scale == "odds") aux$scale <- "log_cumodds" 
    if(aux$scale == "normal") aux$scale <- "inv_normal" 
  } else{
    aux <- NULL
  } 
  return(aux)
}

#' @export
#' @rdname create_params
create_params.flexsurvreg <- function(object, n = 1000, point_estimate = FALSE, ...){
  if (point_estimate == FALSE){
      sim <- flexsurv::normboot.flexsurvreg(object, B = n, raw = TRUE, 
                                            transform = TRUE)
      n_samples <- n
  } else{
      sim <- t(object$res.t[, "est", drop = FALSE])
      n_samples <- 1
  }
  n_pars <- length(object$dlist$pars)
  coefs <- vector(length = n_pars, mode = "list")
  names(coefs) <- c(object$dlist$pars)
  inds <- flexsurvreg_inds(object)
  for (j in seq_along(object$dlist$pars)){
    coefs[[j]] <- sim[, inds[[j]], drop = FALSE]
  }
  return(new_params_surv(dist = object$dlist$name,
                         coefs = coefs,
                         n_samples = n_samples,
                         aux = flexsurvreg_aux(object)))
}

#' @export
#' @rdname create_params
create_params.flexsurvreg_list <- function(object, n = 1000, point_estimate = FALSE, ...){
  return(create_params_list(object, n, point_estimate, 
                          inner_class = "params_surv", new_class = "params_surv_list"))
}

#' @export
#' @rdname create_params
create_params.partsurvfit <- function(object, n = 1000, point_estimate = FALSE, 
                                      bootstrap = TRUE, max_errors = 0, ...){
  if (point_estimate == TRUE & bootstrap == TRUE){
    msg <- paste0("When 'point_estimate' = TRUE and 'bootstrap' = TRUE, the 'point_estimate' ",
                  "argument is ignored and bootstrap replications are generated.")
    warning(msg, call. = FALSE)
  }
  if(bootstrap){
    res <- bootstrap(object, B = n, max_errors = max_errors)
  } else{
    res <- create_params(object$models, n = n, point_estimate = point_estimate)
  }
  return(res)
}

#' @export
# #' @rdname create_params
#' @keywords internal
create_params.joined_flexsurvreg <- function(object, n = 1000, point_estimate = FALSE, ...){
  return(create_params_joined(object, n, point_estimate, "params_surv"))
}

#' @export
# #' @rdname create_params
#' @keywords internal
create_params.joined_flexsurvreg_list <- function(object, n = 1000, point_estimate = FALSE, ...){
  return(create_params_joined(object, n, point_estimate, "params_surv_list"))
}