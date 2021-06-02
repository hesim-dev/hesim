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

# Generic documentation for summary method -------------------------------------
#' Summarize parameter objects
#' 
#' Summarize the coefficients of a parameter object by computing the point estimate, 
#' lower confidence limit, and upper confidence limit for each model term. The point 
#' estimate is the mean of the samples of the coefficients and the lower and 
#' upper confidence limits are determined by the `prob` argument. This is a 
#' convenient way to check whether a parameter object has been specified 
#' correctly and sampling distributions of the coefficients are as expected.
#' @param object An object of the appropriate class.
#' @param prob A numeric scalar in the interval `(0,1)` giving the confidence 
#' interval for coefficients. Default is 0.95 for a 95 percent interval, in which case
#' the lower and upper limits are computed using the 2.5th and 97.5th percentiles.
#' @param ... Additional arguments affecting the summary. Currently unused. 
#' 
#' @return A [`data.table`] containing at least the following columns:
#' \describe{
#' \item{term}{The regression term.}
#' \item{estimate}{The estimated value of the regression term, computed as the mean
#' from its probability distribution.}
#' \item{lower}{The lower limit of the confidence interval for the estimate.}
#' \item{upper}{The upper limit of the confidence interval for the estimate.}
#' }
#' In addition, the following columns may also be present:
#' \describe{
#' \item{parameter}{The name of the parameter of interest. This is relevant
#' for any parametric model in which the underlying probability distribution
#' has multiple parameters. For instance, both [`params_surv`] and [`params_surv_list`]
#' store regression coefficients that are used to model the underlying parameters 
#' of the survival distribution (e.g., shape and scale for a Weibull model). Similarly,
#' there are two parameters (`mean` and `sd`) for [`params_lm`] objects.}
#' \item{model}{The name of the statistical model. This is used for a
#' [`params_surv_list`] object, where each list element represents a separate model.
#' In a state transition model, each model is a unique health state transition and
#' in a partitioned survival model, there is a separate model for each curve.}
#' \item{to}{The health state that is being transitioned to. In [`params_mlogit`]
#'  and [`params_mlogit_list`] objects, there are coefficients for each health
#'  state that can be transitioned to.}
#'  \item{from}{The health state that is being transitions from. This is used
#' for a [`params_mlogit_list`] objects where a different multinomial 
#' logistic regression is used for each state that can be transitioned from.}
#' }
#' 
#' @seealso For examples, see the the underlying parameter object functions: 
#' [params_surv()], [params_surv_list()], [params_mlogit()], and 
#' [params_mlogit_list()].
#' @name summary.params
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
#' distribution. If `"bootstrap"`, then parameters are bootstrapped using [`bootstrap`].
#' If `"none"`, then only point estimates are returned.
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
#' fit <- lm(costs ~ female, data = psm4_exdata$costs$medical)
#' n <- 5
#' params_lm <- create_params(fit, n = n)
#' head(params_lm$coefs)
#' head(params_lm$sigma)
#' 
#' # create_params.flexsurvreg
#' library("flexsurv")
#' fit <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                    data = ovarian, dist = "weibull")
#' n <- 5
#' params_surv_wei <- create_params(fit, n = n)
#' print(params_surv_wei$dist)
#' head(params_surv_wei$coefs)
#' @export
#' @seealso These methods are typically used alongside [create_input_mats()] 
#' to create model objects as a function of input data and a 
#' fitted statistical model. For instance, 
#' [create_PsmCurves()] creates a partitioned survival model,
#' [create_IndivCtstmTrans()] creates an individual continuous time
#' state transition model, 
#' [create_CohortDtstmTrans()] creates a cohort discrete time state transition 
#' model, and
#' [`create_StateVals()`] creates a health state values model. 
#' @rdname create_params
create_params <- function (object, ...) {
  UseMethod("create_params", object)
}

