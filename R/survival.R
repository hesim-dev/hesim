#' Randomly sample parameters from a fitted model
#' 
#' \code{rsample} is a generic function for randomly sampling parameters from
#' a fitted statistical model. Parameters are sampled using the 
#' multivariate normal distribution.
#' @param x A statistical model to randomly sample paraameters from.
#' @param n Number of random observations to draw.
#' @param ... Additional arguments to pass to random sampling functions.
#' @return The form of \code{rsample} depends on the class of its argument. See the 
#' documentation of the particular methods for details of what is produced.
#' @seealso \link{rsample.flexsurvreg}, \link{rsample.flexsurvspline}
#' @name rsample
#' @export
rsample <- function (x, n, ...) {
  UseMethod("rsample", x)
}

rsample_flexsurv <- function(x, n){
  sim <- flexsurv::normboot.flexsurvreg(x, B = n, raw = TRUE)
  n.pars <- length(x$dlist$pars)
  par.samples <- vector(length = n.pars, mode = "list")
  names(par.samples) <- x$dlist$pars
  for (j in seq_along(x$dlist$pars)){
    parname.j <-  x$dlist$pars[j]
    covind.j <- x$mx[[parname.j]]
    if (length(covind.j) > 0){
      ind.j <- c(j, n.pars + covind.j)
    } else{
      ind.j <- j
    }
    par.samples[[j]] <- sim[, ind.j, drop = FALSE]
  }
  return(par.samples)
}

#' Randomly sample parameters from a parametric survival model
#'
#' Draw random samples from a survival model fit using 
#' \code{flexsurv::flexsurvreg}.
#'
#' @param n Number of random observations to draw.
#' @param x Output from \link{flexsurvreg} from the package \link{flexsurv}.
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' @return List where each element contains random samples of the coefficients
#' used to predict a given parameter of the survival distribution. 
#' @name rsample.flexsurvreg
#' @examples 
#' library("flexsurv")
#' fit <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                    data = ovarian, dist = "weibull")
#' n <- 5
#' fit.sample <- rsample(fit, n = n)
#' head(fit.sample)
#' @export
rsample.flexsurvreg <- function(x, n, ...){
  return(rsample_flexsurv(x, n))
}

#' Randomly sample parameters from a spline-based survival model
#'
#' Draw random samples from a survival model fit using 
#' \code{flexsurv::flexsurvspline}.
#'
#' @param n Number of random observations to draw.
#' @param x Output from \link{flexsurvspline} from the package \link{flexsurv}.
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' @return List where each element contains random samples of the coefficients
#' used to predict a given parameter of the survival distribution. 
#' @name rsample.flexsurvspline
#' @examples 
#' library("flexsurv")
#' fit <- flexsurv::flexsurvspline(Surv(recyrs, censrec) ~ group, data = bc, k=1, 
#'                      scale = "hazard")
#' fit.sample <- rsample(fit, n = 2)
#' head(fit.sample)
#' @export
rsample.flexsurvspline <- function(x, n, ...){
  return(rsample_flexsurv(x, n))
}


