# params_surv_list object ------------------------------------------------------
#' Parameters of a list of survival models
#' 
#' Create a list containing the parameters of multiple fitted parametric survival models.
#' @param ... Objects of class [`params_surv`], which can be named.
#' 
#' @return An object of class `params_surv_list`, which is a list containing [`params_surv`]
#' objects.
#' @examples 
#' n <- 5
#' params <- params_surv_list(
#'   # Model 1
#'   mod1 = params_surv(
#'     coefs = list(rate = data.frame(intercept = rlnorm(n, 0, 1))),
#'     dist = "exp"
#'   ),
#'  
#'   # Model 2
#'   mod2 = params_surv(
#'     coefs = list(rate = data.frame(intercept = rlnorm(n, 1, 1))),
#'     dist = "exp"
#'   )
#' )
#' summary(params)
#' params
#' @seealso [create_params()]
#' @export
params_surv_list <- function(...){
  return(check(new_params_list(..., inner_class = "params_surv", 
                               new_class = "params_surv_list")))
}

#' @rdname check
check.params_surv_list <- function(object){
  check_params_list(object)
}

# summary.params_surv_list() ---------------------------------------------------
#' @export
#' @rdname summary.params
summary.params_surv_list <- function(object, prob = 0.95, ...) {
  res <- lapply(object, summary)
  rbindlist(res, idcol = "model")
}

# print.params_surv_list() ----------------------------------------------------------
#' @export
print.params_surv_list <- function(x, ...) {
  
  # Standard output
  cat("A \"params_surv_list\" object:\n\n")
  cat("Summary of coefficient estimates:\n")
  print(summary(x))
  cat("\n")
  cat(paste0("Number of parameter samples: ", x[[1]]$n_samples))
  cat("\n")
  dists <- sapply(x, function(z) z$dist)
  cat("Distributions: ")
  cat("\n")
  print(dists)
  
  # Auxiliary arguments
  if (any(dists %in% c("survspline", "fracpoly", "pwexp"))) {
    cat("Inspect each element of the list to view values for auxilliary parameters.")
  }
  
  # Invisible return
  invisible(x)
}

# create_params methods --------------------------------------------------------
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