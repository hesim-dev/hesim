#' List of \code{formula} objects 
#'
#' Create an object of class "formula_list". The object can be created from either
#'  multiple objects of class \code{\link{formula}} or from another "formula_list" object.
#'
#' @param ... Objects of class \code{\link{formula}}, which can be named.
#' @return Returns an object of class "formula_list".
#' @examples 
#' # Create from "formula" objects
#' formula.list.wei <- formula_list(shape = formula(~ 1), scale = formula(~ x))
#' class(formula.list.wei)
#' 
#' # Create from "formula_list" objects
#' formula.list <- formula_list(exponential = formula_list(rate = formula(~1)),
#'                               weibull = formula.list.wei)
#' 
#' @export
formula_list <- function(...){
  if (inherits(form_object_list(...)[[1]], "formula")){
    return(object_list(..., inner_class = "formula", new_class = "formula_list")) 
  } else {
    return(object_list(..., inner_class = "formula_list", new_class = "formula_list")) 
  } 
}

#' List of \code{lm} objects 
#'
#' Return an object of class "lm_list" multiple objects of class
#' \code{\link{lm}}.
#' @param ... Objects of class \code{\link{lm}}, which can be named.
#' @return Returns an object of class "lm_list".
#' @keywords internal
#' @export
#' @examples 
#'  dat <- part_surv4_simdata$costs$medical
#'  lm.list <- lm_list(fit1 = stats::lm(costs ~ 1, data = dat), 
#'                     fit2 = stats::lm(costs ~ female, data = dat))
#'  class(lm.list)
lm_list <- function(...){
  return(object_list(..., inner_class = "lm", new_class = "lm_list"))
}

#' List of \code{flexsurvreg} objects 
#'
#' Return an object of class "flexsurvreg_list" from multiple objects of class
#' \code{\link{flexsurvreg}}.
#' @param ... Objects of class \code{\link{flexsurvreg}}, which can be named.
#' @return Returns an object of class "flexsurvreg_list".
#' @examples 
#'  library("flexsurv")
#'  fit1 <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist = "weibull")
#'  fit2 <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist = "exp")
#'  flexsurvreg.list <- flexsurvreg_list(wei = fit1, exp = fit2)
#'  class(flexsurvreg.list)
#' @export
flexsurvreg_list <- function(...){
  return(object_list(..., inner_class = "flexsurvreg", new_class = "flexsurvreg_list"))
}

#' Partitioned survival regression object
#'
#' Create a partitioned survival regression object of class "partsurvfit". The object contains a list
#' of fitted survival models fit using either \code{\link{flexsurvreg}} or \code{\link{flexsurvspline}} (i.e.,
#' an object of class \code{\link{flexsurvreg_list}}) and the data frame used to peform the fit of each model. 
#' The same data frame must have been used for each fit.  
#' @param object An object of class \code{\link{flexsurvreg_list}}.
#' @param data The data frame used to fit each survival model in \code{object}.
#' \code{\link{flexsurvreg}}.  
#' @return Returns an object of class "partsurvfit", which is a list containing two elements. 
#' The first element, "models", contains the survival models passed to \code{object}, and the second 
#' element, "data" contains the data frame passed to \code{data}.
#' @examples 
#' library("flexsurv")
#' fit1 <- flexsurv::flexsurvreg(formula = Surv(endpoint1_time, endpoint1_status) ~ age, 
#'                               data = part_surv4_simdata$survival,
#'                               dist = "weibull")
#' fit2 <- flexsurv::flexsurvreg(formula = Surv(endpoint2_time, endpoint2_status) ~ age, 
#'                               data = part_surv4_simdata$survival, 
#'                               dist = "weibull")
#' flexsurvreg.list <- flexsurvreg_list(endpoint1 = fit1, endpoint2 = fit2)
#' fits <- partsurvfit(flexsurvreg.list, data = part_surv4_simdata$survival)
#' class(fits)
#' @export
partsurvfit <- function(object, data){
  if(!inherits(object, "flexsurvreg_list")){
    stop("'Object' must be of class 'flexsurvreg_list'.")
  }
  stopifnot(is.data.frame(data) | is.data.table(data))
  res <- list(models = object, data = data)
  class(res) <- "partsurvfit"
  return(res)
}

#' Join statistical models at specified times
#'
#' Construct classes which are used to join statistical models at specified time points.
#' The first model is used prior to the first time point, the second
#' model is used prior to the second time point, and so on. There should be \code{N + 1}
#' models and \code{N} time points. Can be used to join a single set of models (e.g., for survival analysis or for making predictions from a jointly
#' estimated multi-state model) or to join multiple sets of models (e.g., predicting curves in a 
#' partitioned survival analysis).
#'
#' @param ... Model objects to join. Currently supports objects of class
#' \code{\link{formula}}, \code{\link{formula_list}}, \code{\link{flexsurvreg}}, and
#' \code{\link{flexsurvreg_list}}.
#' @param times Times at which to join models. If joining a single set of models, then a sorted
#' numeric vector of times. If joining multiple sets of models, then a list of sorted numeric vectors,
#' with the length of each list element equal to the number of sets of models. 
#' @return Returns a class of the same name as the function that is a 
#' list with two variables:
#' \describe{
#' \item{models}{List of joined models.}
#' \item{times}{Time points to join models.}
#' }
#' @keywords internal
#' @export
#' @examples 
#' # joined_formula
#' joined.formula <- joined_formula(f1 = formula(~ age), f2 = formula(~ 1),
#'                                  f3 = formula(~age + age2),
#'                                  times = c(2, 3))
#' print(joined.formula)
#'  
#' # joined_formula_list
#' f.list1 <- formula_list(f1 = formula(~ age), f2 = formula(~ 1))
#' f.list2 <- formula_list(f1 = formula(~ age), f2 = formula(~ age2), f3 = formula(~age3))
#' joined.formula.list <- joined_formula_list(model1 = f.list1, model2 = f.list2, 
#'                                            times = list(5, c(2, 4)))
#' class(joined.formula.list)                                             
#'                                             
#' # joined_flexsurvreg
#' library("flexsurv")
#' fit.wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                                 data = ovarian, dist = "weibull")
#' fit.exp <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                                  data = ovarian, dist = "exp")
#' joined.flexsurvreg <- joined_flexsurvreg(fit.wei, fit.exp, times = 5)
#' class(joined.flexsurvreg)
#' @name joined
joined_formula <- function(..., times){
  return(joined_object(..., times = times, inner_class = "formula",
                       new_class = "joined_formula"))
}

#' @rdname  joined
#' @export
joined_formula_list <- function(..., times){
  return(joined_object(..., times = times, inner_class = "formula_list", 
                       new_class = "joined_formula_list", model_list = TRUE))
}

#' @rdname joined
#' @export
joined_flexsurvreg <- function(..., times){
  return(joined_object(..., times = times, inner_class = "flexsurvreg",
                       new_class = "joined_flexsurvreg"))
}

#' @rdname joined
#' @export
joined_flexsurvreg_list <- function(..., times){
  return(joined_object(..., times = times, inner_class = "flexsurvreg_list",
                       new_class = "joined_flexsurvreg_list", model_list = TRUE))
}


