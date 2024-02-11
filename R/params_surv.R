# params_surv() ----------------------------------------------------------------
#' Parameters of a survival model
#' 
#' Create a list containing the parameters of a single fitted parametric or 
#' flexible parametric survival model.
#' @param coefs A list of length equal to the number of parameters in the 
#' survival distribution. Each element of the list is a matrix of samples
#'  of the regression coefficients under sampling uncertainty used to predict
#' a given parameter. All parameters are expressed on the real line (e.g.,
#' after log transformation if they are defined as positive). Each element
#' of the list may also be an object coercible to a matrix such as a 
#' `data.frame` or `data.table`.
#' @param dist Character vector denoting the parametric distribution. See "Details".
#' @param aux Auxiliary arguments used with splines, fractional polynomial,
#' or piecewise exponential models. See "Details". 
#' 
#' @return An object of class `params_surv`, which is a list containing `coefs`,
#' `dist`, and `n_samples`. `n_samples` is equal to the 
#'  number of rows in each element of `coefs`, which must be the same. The `coefs`
#'  element is always converted into a list of matrices. The list may also contain
#'  `aux` if a spline, fractional polynomial, or piecewise exponential model is 
#'  used. 
#' 
#' @details 
#' Survival is modeled as a function of \eqn{L} parameters \eqn{\alpha_l}. 
#' Letting \eqn{F(t)} be the cumulative distribution function, survival at time \eqn{t}
#' is given by
#' \deqn{1 - F(t | \alpha_1(x_{1}), \ldots, \alpha_L(x_{L})).}
#' The parameters are modeled as a function of covariates, \eqn{x_l}, with an 
#' inverse transformation function \eqn{g^{-1}()},
#' \deqn{\alpha_l =  g^{-1}(x_{l}^T \beta_l).}
#' \eqn{g^{-1}()} is typically \eqn{exp()} if a parameter is strictly positive
#' and the identity function if the parameter space is unrestricted.
#' 
#' The types of distributions that can be specified are: 
#' \describe{
#' \item{`exponential` or `exp`}{ Exponential distribution. `coef`
#' must contain the `rate` parameter on the log scale and the same parameterization as in 
#' [`stats::Exponential`].}
#' 
#' \item{`weibull` or `weibull.quiet`}{ Weibull distribution. The first 
#' element of `coef` is the `shape` parameter (on the log scale) and the second
#' element is the `scale` parameter (also on the log scale). The parameterization is
#' that same as in [`stats::Weibull`].}
#' 
#' \item{`weibullPH`}{ Weibull distribution with a proportional hazards 
#' parameterization. The first element of `coef` is the `shape` parameter 
#' (on the log scale) and the second element is the `scale` parameter 
#' (also on the log scale). The parameterization is
#' that same as in [`flexsurv::WeibullPH`].}
#' 
#' \item{`gamma`}{ Gamma distribution. The first 
#' element of `coef` is the `shape` parameter (on the log scale) and the second
#' element is the `rate` parameter (also on the log scale). The parameterization is
#' that same as in [`stats::GammaDist`].}
#' 
#' \item{`lnorm`}{ Lognormal distribution. The first 
#' element of `coef` is the `meanlog` parameter (i.e., the mean of survival on 
#' the log scale) and the second element is the `sdlog` parameter (i.e.,
#'  the standard deviation of survival on the log scale). The parameterization is
#' that same as in [`stats::Lognormal`]. The coefficients predicting the `meanlog` 
#' parameter are untransformed whereas the coefficients predicting the `sdlog` 
#' parameter are defined on the log scale.}
#' 
#' \item{`gompertz`}{ Gompertz distribution. The first 
#' element of `coef` is the `shape` parameter and the second
#' element is the `rate` parameter (on the log scale). The parameterization is
#' that same as in [`flexsurv::Gompertz`].}
#' 
#' \item{`llogis`}{ Log-logistic distribution. The first 
#' element of `coef` is the `shape` parameter (on the log scale) and the second
#' element is the `scale` parameter (also on the log scale). The parameterization is
#' that same as in [`flexsurv::Llogis`].}
#' 
#' \item{`gengamma`}{ Generalized gamma distribution. The first 
#' element of `coef` is the location parameter `mu`, the second
#' element is the scale parameter `sigma` (on the log scale), and the
#' third element is the shape parameter `Q`. The parameterization is
#' that same as in [`flexsurv::GenGamma`].}
#' 
#' \item{`survspline`}{ Survival splines. Each element of `coef` is a parameter of the
#' spline model (i.e. `gamma_0`, `gamma_1`, \eqn{\ldots}) with length equal
#' to the number of knots (including the boundary knots). See below for details on the
#' auxiliary arguments. The parameterization is that same as in [`flexsurv::Survspline`].}
#' 
#' \item{`fracpoly`}{ Fractional polynomials. Each element of `coef` is a parameter of the
#' fractional polynomial model (i.e. `gamma_0`, `gamma_1`, \eqn{\ldots}) with length equal
#' to the number of powers plus 1. See below for details on the auxiliary arguments 
#' (i.e., `powers`).}
#' 
#' \item{`pwexp`}{ Piecewise exponential distribution. Each element of `coef` is 
#' rate parameter for a distinct time interval. The times at which the rates 
#' change should be specified with the auxiliary argument `time` (see below
#' for more details)}. 
#' 
#' \item{`fixed`}{ A fixed survival time. Can be used for "non-random" number 
#' generation. `coef` should contain a single parameter, `est`, of the fixed
#' survival times.}
#' }
#' 
#' Auxiliary arguments for spline models should be specified as a list containing the elements:
#' \describe{
#' \item{`knots`}{A numeric vector of knots.}
#' \item{`scale`}{The survival outcome to be modeled
#' as a spline function. Options are `"log_cumhazard"` for the log cumulative hazard; 
#' `"log_hazard"` for the log hazard rate; `"log_cumodds"` for the log cumulative odds;
#' and `"inv_normal"` for the inverse normal distribution function.}
#' \item{`timescale`}{If `"log"` (the default), then survival is modeled as a spline function
#' of log time; if `"identity"`, then it is modeled as a spline function of time.}
#' }
#' 
#' Auxiliary arguments for fractional polynomial models should be specified as a list containing the elements:
#' \describe{
#' \item{`powers`}{ A vector of the powers of the fractional polynomial with each element
#'  chosen from the following set: -2. -1, -0.5, 0, 0.5, 1, 2, 3.}
#' }
#' 
#' Auxiliary arguments for piecewise exponential models should be specified as 
#' a list containing the element:
#' \describe{
#' \item{`time`}{ A vector equal to the number of rate parameters giving the 
#' times at which the rate changes.}
#' }
#' 
#' Furthermore, when splines (with `scale = "log_hazard"`) or fractional 
#' polynomials are used, numerical methods must be used to compute the cumulative 
#' hazard and for random number generation. The following additional auxiliary arguments
#' can therefore be specified:
#' \describe{
#' \item{`cumhaz_method`}{Numerical method used to compute cumulative hazard 
#' (i.e., to integrate the hazard function). Always used for fractional polynomials
#' but only used for splines if `scale = "log_hazard"`.
#' Options are `"quad"` for adaptive quadrature and `"riemann"` for Riemann sum.}
#'  \item{`random_method`}{Method used to randomly draw from
#'  an arbitrary survival function. Options are `"invcdf"` for the inverse CDF and
#'  `"discrete"` for a discrete time approximation that randomly samples
#'   events from a Bernoulli distribution at discrete times.}
#'  \item{\code{step}}{Step size for computation of cumulative hazard with 
#' numerical integration. Only required when using `"riemann"` to compute the 
#' cumulative hazard or using `"discrete"` for random number generation.}
#' }
#' 
#' @examples 
#' n <- 10
#' params <- params_surv(
#'   coefs = list(
#'     shape = data.frame(
#'       intercept = rnorm(n, .5, .23)
#'     ),
#'     scale = data.frame(
#'       intercept = rnorm(n, 12.39, 1.49),
#'       age = rnorm(n, -.09, .023)
#'    )
#'   ),
#'   dist = "weibull"
#' )
#' summary(params)
#' params
#' @aliases print.params_surv
#' @export
params_surv <- function(coefs, dist, aux = NULL){
  coefs <- coeflist(coefs)
  n_samples <- get_n_samples(coefs)
  check(new_params_surv(coefs, dist, n_samples, aux))
}

new_params_surv <- function(coefs, dist, n_samples, aux = NULL){
  stopifnot(is.vector(dist) & is.character(dist))
  stopifnot(is.numeric(n_samples))
  res <- list(coefs = coefs, dist = dist)
  if (!is.null(aux)) {
    
    # Spline defaults
    if (dist == "survspline") {
      if (is.null(aux$knots)) stop("'knots' must be specified in a spline model.",
                                   call. = FALSE)
      if (is.null(aux$scale)) aux$scale <- "log_cumhazard"
      if (is.null(aux$timescale)) aux$timescale <- "log"
    }
    
    # RNG 
    if (dist %in% c("survspline", "fracpoly")){
      if (is.null(aux$random_method)){
        aux$random_method <- "invcdf"
      } else{
        if (aux$random_method == "sample"){
          warning("'random_method' = 'sample' is deprecated. Use 'discrete' instead.",
                  call. = FALSE)
          aux$random_method <- "discrete"
        }
      }
      
      # Cumulative hazard default
      if (is.null(aux$cumhaz_method)){
        if (dist == "survspline"){
          if (aux$scale == "log_hazard"){
            aux$cumhaz_method <- "quad"
          }
        } else {
          aux$cumhaz_method <- "quad"
        }
      }
    } # End RNG for fractional polynomials and survival splines
    
    # Step size warning
    check_step <- function(aux){
      if (aux$random_method == "discrete" | aux$cumhaz_method == "riemann"){
        if (is.null(aux$step)){
          stop("'step' must be specified.", call. = FALSE)
        }  
      }
    }
    if (dist == "survspline"){
      if (aux$scale == "log_hazard"){
        check_step(aux)
      }
    } 
    if (dist == "fracpoly"){
      check_step(aux)
    }
    
    res[["aux"]] <- aux
    
  } # End if statement for aux
  res[["n_samples"]] <- n_samples
  class(res) <- "params_surv"
  return(res)
}

#' @rdname check
check.params_surv <- function(object, ...){
  # Check coefficients
  if (list_depth(object$coefs) !=1 | length(object$dist) !=1){
    stop("'coefs' must only contain one survival model.", call. = FALSE)
  }
  check(object$coefs)
  
  # Check auxiliary arguments
  # Provide nicer error messages that match.arg()
  check_aux_match <- function(arg, arg_name, choices) {
   if (length(arg) > 1) {
     stop(paste0("The auxiliary argument ", arg_name, " must be of length 1."),
          call. = FALSE)
   } 
   if(!arg %in% choices) {
     stop(paste0("The auxiliary argument ", arg_name, " must be one of ",
                 paste(dQuote(choices), collapse = ", "), "."),
          call. = FALSE)
   }
  }
  
  ## Survival spines
  if (object$dist == "survspline") {
    ### Knots
    if (length(object$coefs) != length(object$aux$knots)) {
      stop(paste0("The number of knots (including boundary knots) in a spline-based ",
                  "survival model must equal the number of parameters."),
           call. = FALSE)
    }
    
    ### Scale
    scale_choices <- c("log_cumhazard", "log_hazard", "log_cumodds", "inv_normal")
    check_aux_match(object$aux$scale, "'scale'", scale_choices)
    
    ### Time scale
    check_aux_match(object$aux$timescale, "'timescale'", c("log", "identity"))
  } 
  
  ## Piecewise exponential
  if (object$dist == "pwexp"){
    if (!is_1d_vector(object$aux$time)){
      stop("`time` must be a 1-dimensional vector", call. = FALSE)
    }
    if (length(object$aux$time) != length(object$coefs)){
      stop("The length of 'time' must equal the length of 'coefs'.", 
           call. = FALSE)
    }
  }
  
  ## Fractional polynomials
  if (object$dist == "fracpoly"){
    if (length(object$coefs) != length(object$aux$powers) + 1) {
      stop(paste0("The number of parameters in a fractional polynomial ",
                  "model must equal the number of powers plus 1."),
           call. = FALSE)
    }
  }
  
  
  # Return
  return(object)
}

# summary.params_surv() --------------------------------------------------------
#' @rdname summary.params
#' @export
summary.params_surv <- function(object, probs = c(.025, .975), ...) {
  
  rbindlist(lapply(object$coef, coef_summary, probs = probs),
            idcol = "parameter")
}

# print.params_surv() ----------------------------------------------------------
#' @export
print.params_surv <- function(x, ...) {
  
  # Standard output
  cat("A \"params_surv\" object \n\n")
  cat("Summary of coefficients:\n")
  print(summary(x))
  cat("\n")
  cat(paste0("Number of parameter samples: ", x$n_samples))
  cat("\n")
  cat(paste0("Distribution: ", x$dist))
  
  # Auxiliary arguments
  if (!is.null(x$aux)) {
    a <- list() 
    if (x$dist == "pwexp") {
      a[["Times:"]] <- x$aux$time
    }
    if (x$dist == "survspline") {
      a[["Knots:"]] <- x$aux$knots
      a[["Scale:"]] <- x$aux$scale
      a[["Time scale:"]] <- x$aux$timescale
    }
    if (x$dist == "fracpoly") {
      a[["Powers:"]] <- x$aux$powers
    }
    for (i in 1:length(a)) {
      cat("\n") 
      cat(names(a)[i], a[[i]])
    }
  } # End printing for auxiliary arguments
  
  # Invisible return
  invisible(x)
}

# create_params.flexsurvreg() --------------------------------------------------
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
create_params.flexsurvreg <- function(object, n = 1000, uncertainty = c("normal", "none"),
                                      ...){
  uncertainty <- deprecate_point_estimate(list(...)$point_estimate, uncertainty,
                                          missing(uncertainty))
  uncertainty <- match.arg(uncertainty)
  if (uncertainty == "normal"){
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
    colnames(coefs[[j]])[1] <- "(Intercept)"
  }
  return(new_params_surv(dist = object$dlist$name,
                         coefs = coefs,
                         n_samples = n_samples,
                         aux = flexsurvreg_aux(object)))
}