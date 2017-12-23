#' Randomly sample parameters from a fitted model
#' 
#' \code{rsample} is a generic function for randomly sampling parameters from
#' a fitted statistical model. Parameters are sampled using the 
#' multivariate normal distribution.
#' @param x A statistical model to randomly sample parameters from.
#' @param ... Additional arguments to pass to random sampling functions.
#' @return The form of \code{rsample} depends on the class of its argument. See the 
#' documentation of the particular methods for details of what is produced.
#' @seealso \link{rsample.flexsurvreg}, \link{rsample.list}
#' @name rsample
#' @export
rsample <- function (x, ...) {
  UseMethod("rsample", x)
}

rsample_flexsurv1 <- function(x, n = NULL, point_estimate = FALSE){
  if (point_estimate == FALSE){
      if (is.null(n)){
        stop("If point_estimate is equal to FALSE, then n must be provided.")
      }
      sim <- flexsurv::normboot.flexsurvreg(x, B = n, raw = TRUE, transform = TRUE)
  } else{
      sim <- t(x$res.t[, "est", drop = FALSE])
  }
  n.pars <- length(x$dlist$pars)
  coefs <- vector(length = n.pars, mode = "list")
  names(coefs) <- c(x$dlist$pars)
  for (j in seq_along(x$dlist$pars)){
    parname.j <-  x$dlist$pars[j]
    covind.j <- x$mx[[parname.j]]
    if (length(covind.j) > 0){
      ind.j <- c(j, n.pars + covind.j)
    } else{
      ind.j <- j
    }
    coefs[[j]] <- sim[, ind.j, drop = FALSE]
  }
  return(list(coefs = coefs, dist = x$dlist$name))
}

#' Randomly sample parameters from a flexsurvreg model
#'
#' Draw random samples from a survival model fit using 
#' \code{flexsurv::flexsurvreg} or \code{flexsurv::flexsurvspline}.
#'
#' @param x Output from \link{flexsurvreg} from the package \link{flexsurv}. 
#' @param n Number of random observations to draw.
#' @param point_estimate If TRUE, then the point estimate is used for prediction
#' and no samples are drawn.
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' @return Returns an object of \link{class} "survlist_pars". A list with
#' the following components is returned:
#' \describe{
#' \item{coefs}{List where each element is a matrix of sampled values of 
#' regression coefficients used to predict a parameter of the probability 
#' distribution.}
#' \item{dist}{Name of the probability distribution.}
#' }
#' @name rsample.flexsurvreg
#' @examples 
#' library("flexsurv")
#' fit <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                    data = ovarian, dist = "weibull")
#' n <- 5
#' fit.sample <- rsample(fit, n = n)
#' head(fit.sample$coefs)
#' fit.sample$dist
#' @export
rsample.flexsurvreg <- function(x, n = NULL, point_estimate = FALSE, ...){
  sim <- rsample_flexsurv1(x, n, point_estimate)
  class(sim) <- "surv_pars"
  return(sim)
}

#' Randomly sample parameters from a list of flexsurvreg models
#'
#' Draw random samples from a list of survival models fit using 
#' \code{flexsurv::flexsurvreg} or \code{flexsurv::flexsurvspline}.
#'
#' @param x A list of models fit using \link{flexsurvreg} from the package
#'  \link{flexsurv}. 
#' @param n Number of random observations to draw.
#' @param point_estimate If TRUE, then the point estimate is used for prediction
#' and no samples are drawn.
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' @return Returns an object of \link{class} "survlist_pars", which is a
#'  lists where each element is a list of class "surv_pars" returned by 
#' \link{rsample.flexsurvreg}. 
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
rsample.list <- function(x, n = NULL, point_estimate = FALSE, ...){
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


#' Obtain design matrix used to fit a model
#' 
#' \code{design_matrix} is a generic function for obtaining the design matrix 
#' used to fit a statistical model.
#' @param x A fitted statistical model.
#' @param ... Additional arguments affecting extracted design matrix.
#' @return The form of \code{design_matrix} depends on the class of its argument. See the 
#' documentation of the particular methods for details of what is produced.
#' @seealso \link{design_matrix.flexsurvreg}, \link{design_matrix.list}
#' @name design_matrix
#' @export
design_matrix <- function (x, ...) {
  UseMethod("design_matrix", x)
}

design_matrix_flexsurv1 <- function(x, indices = NULL){
  X <- x$data$mml
  X <- X[match(names(X), x$dlist$pars)]
  for (i in 1:length(X)){
    if (is.null(X[[i]])){
      X[[i]] <- matrix(1, nrow = nrow(x$data$m))
    }
    if (!is.null(indices)){
      X[[i]] <- X[[i]][indices,, drop = FALSE]
    }
  }
  return(X)
}

#' Obtain design matrix used to fit a flexsurvreg model
#'
#' Obtain design matrix used to fit a survival model using
#' \code{flexsurv::flexsurvreg} or \code{flexsurv::flexsurvspline}.
#'
#' @param x Output from \link{flexsurvreg} from the package \link{flexsurv}. 
#' @param indices Indices from desgin matrices to extract. If NULL, each row 
#' from the design matrix using to fit the \link{flexsurvreg} object 
#' is returned. 
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' @return Returns a list of matrices. Each matrix corresponds to a parameter
#' of the probability distribution.
#' @name design_matrix.flexsurvreg
#' @examples 
#' library("flexsurv")
#' fit <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ age, 
#'                    data = ovarian, dist = "weibull")
#' fit.X <- design_matrix(fit)
#' names(fit.X)
#' head(fit.X$scale)
#' head(fit.X$shape)
#' @export
design_matrix.flexsurvreg <- function(x, indices = NULL, ...){
  return(design_matrix_flexsurv1(x, indices))
}

#' Obtain design matrix used to fit a list of flexsurvreg models
#'
#' Obtain design matrix used to fit a survival model using
#' \code{flexsurv::flexsurvreg} or \code{flexsurv::flexsurvspline}.
#'
#' @param x Output from \link{flexsurvreg} from the package \link{flexsurv}. 
#' @param indices Indices from desgin matrices to extract. If NULL, each row 
#' from the design matrix using to fit the \link{flexsurvreg} object 
#' is returned. 
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' @return Returns a list of lists matrices. Each list is equivalent to the
#' output produced by \link{design_matrix.flexsurvreg}.
#' @name design_matrix.list
#' @examples 
#' library("flexsurv")
#' fit1 <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                    data = ovarian, dist = "exponential")
#' fit2 <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                    data = ovarian, dist = "weibull")                   
#' fit.X <- design_matrix(list(fit1, fit2))
#' names(fit.X)
#' head(fit.X[[1]])
#' @export
design_matrix.list <- function(x, indices = NULL, ...){
  ret <- vector(length(x), mode = "list")
  for (i in 1:length(x)){
    if (!inherits(x[[i]], "flexsurvreg")){
      stop("Each element of x must be of class flexsurvreg")
    }
    ret[[i]] <- design_matrix_flexsurv1(x[[i]], indices)
  } 
  return(ret)
}

#' Survival Simulation
#' 
#' An R6 class for simulating health-economic models based on overall survival.
#' @format \code{\link{R6Class}} object.
#' @field dist_name Name of the probability distribution used to model survival.
#' @field  surv_coefs An object of class "surv_pars" representing the posterior
#' distribution of coefficients from a survival model. For example, the output of
#'  \code{\link{rsample}}.
#' @field surv_X Matrix of covariates used to predict survival using the 
#' coefficients from \code{surv_coefs}. 
#' 
#' @section Members:
#'  \describe{
#'   \item{\code{new()}}{Create a new object of class \code{Survival}}
#'   \item{\code{quantiles(q)}}{Generate the quantiles \code{q} of the survival 
#'   distribution. \code{q} is a numeric vector.} 
#'   \item{\code{surv(t)}}{Generate a survival curve give a numeric vector
#'   of times \code{t}.}
#'   \item{\code{cumhazard(t)}}{Generate cumulative hazards give a numeric vector
#'   of times \code{t}.}
#'   \item{\code{hazard(t)}}{Generate hazards give a numeric vector
#'   of times \code{t}.}
#'   \item{\code{rmst(t, r = 0)}}{Calculate the (discounted) restricted mean 
#'   survival time. Calculated for a vector of times \code{t} given a discount
#'   rate \code{r}. Durvival times are not discounted by default.
#' }
#' 
#' @examples 
#' @export
Survival <- R6::R6Class("Survival",
  public = list(
    dist_name = NULL,
    surv_coefs = NULL,
    surv_X = NULL,
    
    initialize = function(dist_name = NULL, surv_coefs = NULL, surv_X = NULL) {
      self$dist_name <- dist_name
      self$surv_coefs <- surv_coefs
      self$surv_X <- surv_X
    },
    
    quantiles = function(q){
      R_survival_summary(self$dist_name, self$surv_coefs, 
                         self$surv_X, 0,
                         q, type = "quantiles")
    },
    
    surv = function(t){
      R_survival_summary(self$dist_name, self$surv_coefs, 
                         self$surv_X, 0, 
                         t, type = "survival")
    },
    
    cumhazard = function(t){
      R_survival_summary(self$dist_name, self$surv_coefs, 
                         self$surv_X, 0,
                         t, type = "cumhazard")
    },
    
    hazard = function(t){
      R_survival_summary(self$dist_name, self$surv_coefs, 
                        self$surv_X, 0,
                         t, type = "hazard")
    },
    
    rmst = function(t, r = 0){
      R_survival_summary(self$dist_name, self$surv_coefs, 
                         self$surv_X, r,
                         t, type = "rmst")
    }
  )
)

R_survival_summary <- function(dist_name, surv_coefs, surv_X, r, 
                               x, type){
  ret <- C_survival_summary(dist_name, surv_coefs, surv_X, r, 
                            x, type)
  ret <- data.table::data.table(ret)
  if (type == "quantiles"){
      setnames(ret, colnames(ret), c("sim", "id", "q", "value"))
  } else{
      setnames(ret, colnames(ret), c("sim", "id", "t", "value"))
  }
  sim = id = NULL # for no visible binding note in CRAN check
  ret[, sim := sim + 1]
  ret[, id := id + 1]
  ret[]
}
