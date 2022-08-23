# params_mlogit_list() ---------------------------------------------------------
#' Parameters of a list of multinomial logit models
#'
#' Create a list containing the parameters of multiple fitted multinomial logit models.
#' Can be used to parameterize state transitions in a discrete time transition model
#' by passing to the `params` field of a [`CohortDtstmTrans`] object.
#' @param ... Objects of class [`params_mlogit`], which can be named.
#'
#' @return An object of class `params_mlogit_list`, which is a list containing
#' [`params_mlogit`] objects.
#' @examples
#' # Consider a sick-sicker model
#'
#' params <- params_mlogit_list(
#'   ## Transitions from sick state (sick -> sicker, sick -> death)
#'   sick = params_mlogit(
#'     coefs = list(
#'       sicker = data.frame(
#'         intercept = c(-0.33, -.2),
#'         treat = c(log(.75), log(.8))
#'       ),
#'       death = data.frame(
#'         intercept = c(-1, -1.2),
#'         treat = c(log(.6), log(.65))
#'       )
#'     )
#'   ),
#'
#'   ## Transitions from sicker state (sicker -> death)
#'   sicker = params_mlogit(
#'     coefs = list(
#'       death = data.frame(
#'         intercept = c(-1.5, -1.4),
#'         treat = c(log(.5), log(.55))
#'       )
#'     )
#'   )
#' )
#' summary(params)
#' params
#'
#' @seealso [summary.params_mlogit_list()], [params_mlogit()], [`CohortDtstmTrans`]
#' @export
params_mlogit_list <- function(...) {
  p <- new_params_list(...,
    inner_class = "params_mlogit",
    new_class = "params_mlogit_list"
  )
  check_params_list(p)
}

# summary.params_surv_list() ---------------------------------------------------
#' @rdname summary.params
#' @export
summary.params_mlogit_list <- function(object, probs = c(.025, .975), ...) {
  summary_params_list(object, probs, idcol = "from", ...)
}

# print.params_mlogit_list() ---------------------------------------------------
#' @export
print.params_mlogit_list <- function(x, ...) {
  cat("A \"params_mlogit_list\" object\n\n")
  cat("Summary of coefficients:\n")
  print(summary(x))
  cat("\n")
  cat(paste0("Number of parameter samples: ", x[[1]]$n_samples))
  cat("\n")
  cat(paste0("Number of starting (non-absorbing) states: ", length(x)))
  cat("\n")
  cat(
    "Number of transitions by starting state:",
    sapply(x, function(z) dim(z$coef)[3])
  )

  invisible(x)
}

# create_params.multinom_list() ------------------------------------------------
#' @export
#' @rdname create_params
create_params.multinom_list <- function(object, n = 1000, uncertainty = c("normal", "none"), ...) {
  return(create_params_list(object,
    n = n, uncertainty = uncertainty,
    inner_class = "params_mlogit", new_class = "params_mlogit_list",
    ...
  ))
}
