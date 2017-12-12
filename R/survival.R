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
#' @seealso \link{rsample.flexsurvreg}, \link{rsample.list}
#' @name rsample
#' @export
rsample <- function (x, n, ...) {
  UseMethod("rsample", x)
}

rsample_flexsurv1 <- function(x, n){
  sim <- flexsurv::normboot.flexsurvreg(x, B = n, raw = TRUE)
  n.pars <- length(x$dlist$pars)
  par.samples <- vector(length = n.pars + 1, mode = "list")
  names(par.samples) <- c(x$dlist$pars, "dist")
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
  par.samples[["dist"]] <- x$dlist$name
  return(par.samples)
}

#' Randomly sample parameters from a flexsurvreg model
#'
#' Draw random samples from a survival model fit using 
#' \code{flexsurv::flexsurvreg} or \code{flexsurv::flexsurvspline}.
#'
#' @param n Number of random observations to draw.
#' @param x Output from \link{flexsurvreg} from the package \link{flexsurv}. 
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' @return Returns an object of \link{class} "surv_pars", which is a list 
#' containing (i) random samples of the regression coefficients used to predict
#' each of the model parameters and (ii) the probability distribution. 
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
  sim <- rsample_flexsurv1(x, n)
  class(sim) <- "surv_pars"
  return(sim)
}

#' Randomly sample parameters from a list of flexsurvreg models
#'
#' Draw random samples from a list of survival models fit using 
#' \code{flexsurv::flexsurvreg} or \code{flexsurv::flexsurvspline}.
#'
#' @param n Number of random observations to draw.
#' @param x A list of models fit using \link{flexsurvreg} from the package
#'  \link{flexsurv}. 
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' @return Returns an object of \link{class} "survlist_pars", which is a list 
#' of lists. Each list contains (i) random samples of the regression coefficients 
#' used to predict each of the model parameters and (ii) the
#'  probability distribution.
#' @name rsample.list
#' @examples 
#' library("flexsurv")
#' fit1 <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ age, 
#'                     data = ovarian, dist = "exponential")
#' fit2 <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                    data = ovarian, dist = "gamma")
#' fits <- list(f1 = fit1, f2 = fit2)
#' fits.sample <- rsample(fits, n = 2)
#' head(fits.sample)
#' @export
rsample.list <- function(x, n, ...){
  ret <- vector(length(x), mode = "list")
  for (i in 1:length(x)){
    if (!inherits(x[[i]], "flexsurvreg")){
      stop("Each element of x must be of class flexsurvreg")
    }
    ret[[i]] <- rsample_flexsurv1(x[[i]], n)
  } 
  class(ret) <- "survlist_pars"
  return(ret)
}

