# Deprecated functions ---------------------------------------------------------
#' Individualized cost-effectiveness analysis
#'
#' These functions are deprecated, use [cea()] and [cea_pw()] instead. 
#' @param x An object of simulation output characterizing the probability distribution
#' of clinical effectiveness and costs.?ic
#' @param ... Further arguments passed to or from other methods. 
#' @export
#' @rdname icea
icea <- function(x, ...) {
  .Deprecated("cea")
  UseMethod("cea")
}

#' @export
#' @rdname icea
icea_pw <- function(x, ...) {
  .Deprecated("cea_pw")
  UseMethod("cea_pw")
}

# Deprecated arguments ---------------------------------------------------------
deprecate_point_estimate <- function(old, new, is_new_missing) {
  if (!is.null(old)) { 
    warning("'point_estimate' is deprecated; use 'uncertainty' instead.",
            call. = FALSE)
  }
  if (!is.null(old) && (old == TRUE & is_new_missing == TRUE)) {
    return("none")
  } else{
    return(new)
  }
}

deprecate_bootstrap <- function(old, new, is_new_missing) {
  if (!is.null(old)) { 
    warning("'bootstrap' is deprecated; use 'uncertainty' instead.",
            call. = FALSE)
  }
  if (!is.null(old) && (old == TRUE & is_new_missing == TRUE)) {
    return("bootstrap")
  } else{
    return(new)
  }
}