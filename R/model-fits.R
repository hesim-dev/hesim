#' List of `formula` objects 
#'
#' Combine [formula][stats::formula] or [formula_list] object into a 
#' [formula_list] object.
#'
#' @param ... Objects of class [formula][stats::formula], which can be named.
#' @return An object of class `formula_list`.
#' @keywords internal
#' @examples 
#' # Create from "formula" objects
#' flist_wei <- formula_list(shape = formula(~ 1), scale = formula(~ x))
#' class(flist_wei)
#' 
#' # Create from "formula_list" objects
#' flist <- formula_list(exponential = formula_list(rate = formula(~1)),
#'                               weibull = flist_wei)
#' 
#' @export
formula_list <- function(...){
  if (inherits(create_object_list(...)[[1]], "formula")){
    return(object_list(..., inner_class = "formula", new_class = "formula_list")) 
  } else {
    return(object_list(..., inner_class = "formula_list", new_class = "formula_list")) 
  } 
}

#' List of `lm` objects 
#'
#' Combine [`lm`][stats::lm] objects into a list
#' @param ... Objects of class [`lm`][stats::lm], which can be named.
#' @return Returns an object of class `lm_list`.
#' @keywords internal
#' @export
#' @examples 
#'  dat <- psm4_exdata$costs$medical
#'  lm_fits <- lm_list(fit1 = stats::lm(costs ~ 1, data = dat), 
#'                     fit2 = stats::lm(costs ~ female, data = dat))
#'  class(lm_fits)
lm_list <- function(...){
  return(object_list(..., inner_class = "lm", new_class = "lm_list"))
}

#' List of `flexsurvreg` objects 
#'
#' Combine [`flexsurvreg`][flexsurv::flexsurvreg] objects into a list.
#' @param ... Objects of class [`flexsurvreg`][flexsurv::flexsurvreg], which can be named.
#' @return An object of class `flexsurvreg_list`.
#' @examples 
#'  library("flexsurv")
#'  fit1 <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist = "weibull")
#'  fit2 <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist = "exp")
#'  fsreg_list <- flexsurvreg_list(wei = fit1, exp = fit2)
#'  class(fsreg_list)
#' @export
flexsurvreg_list <- function(...){
  return(object_list(..., inner_class = "flexsurvreg", new_class = "flexsurvreg_list"))
}

#' List of `multinom` objects 
#'
#' Combine `multinom` objects into a list.
#' @param ... Objects of class [`multinom`][nnet::multinom], which can be named.
#' @return  An object of class `multinom_list`.
#' @examples 
#'  library("nnet")
#'  library("data.table")
#'  trans_data <- data.table(multinom3_exdata$transitions)
#'  dat_healthy <- trans_data[state_from == "Healthy"]
#'  fit_healthy <- multinom(state_to ~ strategy_name + female + age_cat + year_cat, 
#'                           data = dat_healthy)
#'  dat_sick <- trans_data[state_from == "Sick"]
#'  dat_sick$state_to <- droplevels(dat_sick$state_to)
#'  fit_sick <- multinom(state_to ~ strategy_name + female + age_cat + year_cat, 
#'                       data = dat_sick)
#'  fits <- multinom_list(healthy = fit_healthy, sick = fit_sick)
#'  class(fits)
#' @export
multinom_list <- function(...){
  return(object_list(..., inner_class = "multinom", new_class = "multinom_list"))
}

#' Partitioned survival regression object
#'
#' Create a partitioned survival regression object of class `partsurvfit`. The object contains a list
#' of fitted survival models fit using either [`flexsurv::flexsurvreg`]
#' or [`flexsurv::flexsurvspline`] (i.e.,
#' an object of class \code{\link{flexsurvreg_list}}) and the data frame used to perform the fit of each model. 
#' The same data frame must have been used for each fit.  
#' @param object An object of class [`flexsurv::flexsurvreg_list`].
#' @param data The data frame used to fit each survival model in `object`.
#' @return Returns an object of class `partsurvfit`, which is a list containing two elements. 
#' The first element, "models", contains the survival models passed to `object`, and the second 
#' element, "data" contains the data frame passed to `data`.
#' @examples 
#' library("flexsurv")
#' fit1 <- flexsurv::flexsurvreg(formula = Surv(endpoint1_time, endpoint1_status) ~ age, 
#'                               data = psm4_exdata$survival,
#'                               dist = "weibull")
#' fit2 <- flexsurv::flexsurvreg(formula = Surv(endpoint2_time, endpoint2_status) ~ age, 
#'                               data = psm4_exdata$survival, 
#'                               dist = "weibull")
#' fsreg_list <- flexsurvreg_list(endpoint1 = fit1, endpoint2 = fit2)
#' fits <- partsurvfit(fsreg_list, data = psm4_exdata$survival)
#' class(fits)
#' @export
#' @keywords internal
partsurvfit <- function(object, data){
  if(!inherits(object, "flexsurvreg_list")){
    stop("'Object' must be of class 'flexsurvreg_list'.")
  }
  stopifnot(is.data.frame(data) | is.data.table(data) | is.null(data))
  res <- list(models = object, data = data)
  class(res) <- "partsurvfit"
  return(res)
}