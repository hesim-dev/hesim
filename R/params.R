new_params_list <- function(..., inner_class, new_class){
  return(object_list(..., inner_class = inner_class,
                     new_class = new_class))
}

check_params_list <- function(x){
  inner.class <- class(x[[1]])
  n.samples <- sapply(x, function(y) y$n_samples)
  if(!all(n.samples == n.samples[1])){
    msg <- paste0("The number of samples in each '", inner.class, "' object must be the same.")
    stop(msg, call. = FALSE)
  }
  return(x)
}

new_params_joined <- function(..., times, inner_class){
  objects <- form_object_list(...)
  res <- list(models = objects, times = times)
  class(res) <- paste0("params_joined", sub("params", "", inner_class))
  return(res)
}

check_params_joined <- function(x, inner_class, model_list){
  check_joined_object(x, inner_class = inner_class, model_list = model_list)
  return(x)
}

#' Parameters of a linear model
#' 
#' Create a list containing the parmeters of a fitted linear regression model.
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
#' params.lm <- params_lm(coefs = MASS::mvrnorm(n, mu = c(.5,.6),
#'                                             Sigma = matrix(c(.05, .01, .01, .05), nrow = 2)),
#'                       sigma <- rgamma(n, shape = .5, rate = 4))
#' print(params.lm)
#'
#' @export
params_lm <- function(coefs, sigma = NULL){
  stopifnot(is.matrix(coefs) | is.vector(coefs))
  if (is.vector(coefs)) coefs <- matrix(coefs, ncol = 1)
  if(is.null(sigma)) sigma <- rep(1, nrow(coefs))
  n.samples <- nrow(coefs)
  check(new_params_lm(coefs, sigma, n.samples))
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
#' Create a list containing the parmeters of a list of fitted linear regression models.
#' @param ... Objects of class \code{\link{params_lm}}, which can be named.
#' 
#' @return An object of class "params_lm_list", which is a list containing \code{\link{params_lm}}
#' objects. 
#' @export
#' @keywords internal
#' @examples 
#' n <- 2
#' coefs1 <- MASS::mvrnorm(n, mu = c(.5,.6),
#'                        Sigma = matrix(c(.05, .01, .01, .05), nrow = 2))
#' sigma1 <- rgamma(n, shape = .5, rate = 4)
#' params.lm1 <- params_lm(coefs = coefs1, sigma = sigma1)
#' 
#' coefs2 <- MASS::mvrnorm(n, mu = c(.2,.9),
#'                        Sigma = matrix(c(.08, .02, .02, .08), nrow = 2))
#' sigma2 <- rgamma(n, shape = .9, rate = 4)
#' params.lm2 <- params_lm(coefs = coefs2, sigma = sigma2)
#' 
#' params.lm.list <- params_lm_list(params.lm1, params.lm2)
#' print(params.lm.list)
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
#' @param aux Auxillary arguments used with splines. See "Details". 
#' 
#' @return An object of class "params_surv", which is a list containing \code{coefs},
#' \code{dist}, and \code{n_samples}. \code{n_samples} is equal to the number of rows
#' in each element of \code{coefs}, which must be the same. The list may also contain \code{aux}, if
#' a spline model is fit. 
#' 
#' @details 
#' The types of distributions that can be specified are: 
#' \itemize{
#' \item{\code{exponential} or \code{exp}}{ Exponential distribution}
#' \item{\code{weibull} or \code{weibull.quiet}}{ Weibull distribution}
#' \item{\code{gamma}}{ Gamma distribution}
#' \item{\code{lnorm}}{ Lognormal distribution}
#' \item{\code{gompertz}}{ Gompertz distribution}
#' \item{\code{llogis}}{ Log-logistic distribution}
#' \item{\code{gengamma}}{ Generalized gamma distribution}
#' \item{\code{survspline}}{ Survival splines}
#' }
#' 
#' Auxillary arguments for spline models should be specified as a list containing the elements:
#' \describe{
#' \item{\code{knots}}{ A numeric vector of knots.}
#' \item{\code{scale}}{ A character vector of length 1 denoting the survival outcome to be modeled
#' as a spline function. Options are "log_cumhazard", for the log cummulative hazard; 
#' "log_hazard" for the log hazard rate; "log_cumodds" for the log cummulative odds;
#' and "inv_normal" for the inverse normal distribution function.}
#' \item{\code{timescale}}{If "log" (the default), then survival is modeled as a spline function
#' of log time; if "identity", then it is modeled as a spline function of time.}
#' }
#' 
#' @examples 
#' library("flexsurv")
#' fit <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist = "weibull")
#' params.surv <- params_surv(coefs = list(shape = fit$res.t["shape", "est", drop = FALSE],
#'                                         scale = fit$res.t["scale", "est", drop = FALSE]),
#'                            dist = fit$dlist$name)
#' print(params.surv)
#' @export
params_surv <- function(coefs, dist, aux = NULL){
  stopifnot(is.list(coefs))
  if (!is.matrix(coefs[[1]])){
    stop("'coefs' must be a list of matrices.", call. = FALSE)
  }
  n.samples <- nrow(coefs[[1]])
  check(new_params_surv(coefs, dist, n.samples, aux))
}

new_params_surv <- function(coefs, dist, n_samples, aux = NULL){
  stopifnot(is.vector(dist) & is.character(dist))
  stopifnot(is.numeric(n_samples))
  res <- list(coefs = coefs, dist = dist)
  if (!is.null(aux)) res[["aux"]] <- aux
  res[["n_samples"]] <- n_samples
  class(res) <- "params_surv"
  return(res)
}

#' @rdname check
check.params_surv <- function(object){
  if (list_depth(object$coefs) !=1 | length(object$dist) !=1){
    stop("'coefs' must only contain one survival model.", call. = FALSE)
  }
  matrix.bool <- unlist(lapply(object$coefs, is.matrix))
  if(sum(!matrix.bool) > 0){
    stop("'coefs' must be a list of matrices.",
         call. = FALSE)
  }
  coefs.nrows <- unlist(lapply(object$coefs, nrow))
  if(!all(object$n_samples == coefs.nrows)){
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
#' fit.exp <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                                  data = ovarian, dist = "exp")
#' fit.wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                                  data = ovarian, dist = "weibull")
#' params.surv.exp <- form_params(fit.exp, n = 2)
#' params.surv.wei <- form_params(fit.wei, n = 2)
#' params.joined.surv <- params_joined_surv(exp = params.surv.exp,
#'                                          wei = params.surv.wei,
#'                                          times = 3)
#' print(params.joined.surv)
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
#' fit.wei <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist = "weibull")
#' params.surv.wei <- form_params(fit.wei, n = 2)
#' 
#' fit.exp <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist = "exp")
#' params.surv.exp <- form_params(fit.exp, n = 2)
#' 
#' params.surv.list <- params_surv_list(wei = params.surv.wei, exp = params.surv.exp)
#' print(params.surv.list)
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
#' at specified time points.See \code{\link{joined}} for more details.
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
#' fit.exp <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                                  data = ovarian, dist = "exp")
#' fit.wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                                  data = ovarian, dist = "weibull")
#' fit.lnorm <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                                   data = ovarian, dist = "lognormal")
#'
#' params.surv.exp <- form_params(fit.exp, n = 2)
#' params.surv.wei <- form_params(fit.wei, n = 2)
#' params.surv.lnorm <- form_params(fit.lnorm, n = 2)
#'
#' params.surv.list1 <- params_surv_list(params.surv.exp, params.surv.wei)
#' params.surv.list2 <- params_surv_list(params.surv.exp, params.surv.lnorm)
#' params.joined.surv.list <- params_joined_surv_list(model1 = params.surv.list1,
#'                                                    model2 = params.surv.list2,
#'                                                    times = list(3, 5))
#' print(params.joined.surv.list)
#' @export
params_joined_surv_list <- function(..., times){
  return(check(new_params_joined(..., times = times, inner_class = "params_surv_list"),
               inner_class = "params_surv_list"))
}

#' @rdname check
check.params_joined_surv_list <- function(object, inner_class){
  check_params_joined(object, inner_class = inner_class, model_list = TRUE)
}

form_params_list <- function(object, n, point_estimate, inner_class, new_class){
  n.objects <- length(object)
  params.list <- vector(mode = "list", length = n.objects)
  names(params.list) <- names(object)
  for (i in 1:n.objects){
    params.list[[i]] <- form_params(object[[i]], n, point_estimate)
  }
  return(new_params_list(params.list, inner_class = inner_class,
                         new_class = new_class))
}

form_params_joined <- function(object, n, point_estimate, inner_class){
  n.models <- length(object$models)
  models <- vector(mode = "list", length = n.models)
  names(models) <- names(object$models)
  for (i in 1:n.models){
    models[[i]] <- form_params(object$models[[i]], n, point_estimate)
  }
  return(new_params_joined(models, times = object$times, inner_class = inner_class))
}

#' Parameters from a fitted model
#' 
#' \code{form_params} is a generic function for constructing parameters from
#' a fitted statistical model. A sample from the posterior distribution of 
#' parameters is returned.
#' @param object A statistical model to randomly sample parameters from.  
#' @param n Number of random observations to draw. 
#' @param point_estimate If TRUE, then the point estimates are returned and
#' and no samples are drawn.
#' @param bootstrap If \code{bootstrap} is FALSE or not specified, then \code{n} parameter sets are 
#' drawn by sampling from a multivariate normal distribution. If \code{bootstrap} is TRUE, then 
#' parameters are bootstrapped using \code{\link{bootstrap}}. 
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' @return An object prefixed by \code{params_}. Mapping between \code{form_params} 
#' and the classes of the returned objects are: 
#' \itemize{
#' \item{\code{form_params.lm} ->}{ \code{params_lm}}
#' \item{\code{form_params.flexsurvreg} ->}{ \code{params_surv}}
#' \item{\code{form_params.flexsurvreg_list} ->}{ \code{params_surv_list}}
#' \item{\code{form_params.partsurvfit} ->}{ \code{params_surv_list}}
#' }
#' @name form_params
#' @examples 
#' # form_params.lm
#' fit <- stats::lm(costs ~ female, data = part_surv4_simdata$costs$medical)
#' n <- 5
#' params.lm <- form_params(fit, n = n)
#' head(params.lm$coefs)
#' head(params.lm$sigma)
#' 
#' # form_params.flexsurvreg
#' library("flexsurv")
#' fit <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                     data = ovarian, dist = "weibull")
#' n <- 5
#' params.surv <- form_params(fit, n = n)
#' print(params.surv$dist)
#' head(params.surv$coefs)
#' @export
#' @rdname form_params
form_params <- function (object, n = 1000, point_estimate = FALSE, ...) {
  UseMethod("form_params", object)
}

#' @export
#' @rdname form_params
form_params.lm <- function(object, n = 1000, point_estimate = FALSE, ...){
  if (point_estimate == FALSE){
      coefs.sim <- MASS::mvrnorm(n, stats::coef(object), stats::vcov(object))
      if(is.vector(coefs.sim)) coefs.sim <- matrix(coefs.sim, nrow = 1)
      return(new_params_lm(coefs = coefs.sim,
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
# #' @rdname form_params
#' @keywords internal
form_params.lm_list <- function(object, n = 1000, point_estimate = FALSE, ...){
  return(form_params_list(object, n, point_estimate, 
                          inner_class = "params_lm", new_class = "params_lm_list"))
}

flexsurvreg_inds <- function(object){
  n.pars <- length(object$dlist$pars)
  inds <- vector(mode = "list", length = n.pars)
  for (j in 1:n.pars){
    parname.j <-  object$dlist$pars[j]
    covind.j <- object$mx[[parname.j]]
    if (length(covind.j) > 0){
      inds[[j]] <- c(j, n.pars + covind.j)
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
#' @rdname form_params
form_params.flexsurvreg <- function(object, n = 1000, point_estimate = FALSE, ...){
  if (point_estimate == FALSE){
      sim <- flexsurv::normboot.flexsurvreg(object, B = n, raw = TRUE, 
                                            transform = TRUE)
      n.samples <- n
  } else{
      sim <- t(object$res.t[, "est", drop = FALSE])
      n.samples <- 1
  }
  n.pars <- length(object$dlist$pars)
  coefs <- vector(length = n.pars, mode = "list")
  names(coefs) <- c(object$dlist$pars)
  inds <- flexsurvreg_inds(object)
  for (j in seq_along(object$dlist$pars)){
    coefs[[j]] <- sim[, inds[[j]], drop = FALSE]
  }
  return(new_params_surv(dist = object$dlist$name,
                         coefs = coefs,
                         n_samples = n.samples,
                         aux = flexsurvreg_aux(object)))
}

#' @export
#' @rdname form_params
form_params.flexsurvreg_list <- function(object, n = 1000, point_estimate = FALSE, ...){
  return(form_params_list(object, n, point_estimate, 
                          inner_class = "params_surv", new_class = "params_surv_list"))
}

#' @export
#' @rdname form_params
form_params.partsurvfit <- function(object, n = 1000, point_estimate = FALSE, 
                                      bootstrap = TRUE, ...){
  if (point_estimate == TRUE & bootstrap == TRUE){
    msg <- paste0("When 'point_estimate' = TRUE and 'bootstrap' = TRUE, the 'point_estimate' ",
                  "argument is ignored and bootstrap replications are generated.")
    warning(msg, call. = FALSE)
  }
  if(bootstrap){
    res <- bootstrap(object, B = n)
  } else{
    res <- form_params(object$models, n = n, point_estimate = point_estimate)
  }
  return(res)
}

#' @export
# #' @rdname form_params
#' @keywords internal
form_params.joined_flexsurvreg <- function(object, n = 1000, point_estimate = FALSE, ...){
  return(form_params_joined(object, n, point_estimate, "params_surv"))
}

#' @export
# #' @rdname form_params
#' @keywords internal
form_params.joined_flexsurvreg_list <- function(object, n = 1000, point_estimate = FALSE, ...){
  return(form_params_joined(object, n, point_estimate, "params_surv_list"))
}