#' Random generation for categorical distribution
#'
#' Draw random samples from a categorical distribution given a matrix of probabilities.
#'  \code{rcat} is vetorized and written in C++ for speed.  
#'
#' @param prob A matrix of probabilities where rows correspond to observations 
#' and columns correspond to categories.
#' @name rcat
#' @examples
#' p <- c(.2, .5, .3)
#' n <- 1000000
#' pmat <- matrix(rep(p, n), nrow = n, ncol = length(p), byrow = TRUE)
#'
#' # rcat
#' set.seed(100)
#' ptm <- proc.time()
#' samp1 <- rcat(pmat)
#' proc.time() - ptm
#' prop.table(table(samp1))
#'
#' # rmultinom from base R 
#' set.seed(100)
#' ptm <- proc.time()
#' samp2 <- t(apply(pmat, 1, rmultinom, n = 1, size = 1))
#' samp2 <- apply(samp2, 1, function(x) which(x == 1))
#' proc.time() - ptm
#' prop.table(table(samp2))
#' @export
rcat <- function(prob){
  if(!is.matrix(prob)){
    stop("prob must be a matrix")
  }
  return(rcatC(prob) + 1)
}


#' Random generation for piecewise exponential distribution
#'
#' Draw random samples from an exponential distribution with piecewise rates.
#'  \code{rpwexp} is vetorized and written in C++ for speed.  
#'
#' @param rate A matrix of rates where rows correspond to observations 
#' and columns correspond to rates during specified time intervals. 
#' @param time A vector equal to the number of columns in \code{rate} giving the
#' times at which the rate changes
#' @name rpwexp
#' @examples
#' rate <- c(.6, 1.2, 1.3)
#' n <- 100000
#' ratemat <- matrix(rep(rate, n/2), nrow = n, 
#'                   ncol = 3, byrow = TRUE)
#' t <- c(0, 10, 15) 
#' ptm <- proc.time()
#' samp <- rpwexp(ratemat, t)
#' proc.time() - ptm
#' summary(samp)
#' @export
rpwexp <- function(rate, time){
  return(rpwexpC(rate, time))
}

#' Extract parameters for simulating survival models.
#'
#' Extract parameters for simulating survival or multi-state models fit using \code{flexsurvreg}
#' from the \code{flexsurv} package.
#'
#' @param x For \code{rsurv_prep}, a model fit with \code{flexsurvreg}; for \code{rmsm_prep}, a list of
#' models fit with \code{flexsurvreg}.
#'
#' @return List with the following parameters:
#' \item{loc_beta}{Regression coefficients for location parameter.}
#' \item{dist}{Distribution used for model fit.}
#' \item{par}{1st ancillary parameter in survvial model.}
#'
#' @name rsurv_prep
#' @export
#' @keywords internal
rmsm_prep <- function(x){
  loc_beta <- anc1 <- list()
  dist <- rep(NA, length(x))
  for (i in 1:length(x)){
    param <- rsurv_prep(x[[i]])
    loc_beta[[i]] <- param$location
    anc1[[i]] <- param$anc1
    dist[i] <- param$dist
  }
  loc_beta <- list_to_array(loc_beta)
  anc1 <- unlist(anc1)
  dist <- unlist(dist)
  return(list(loc_beta = loc_beta, anc1 = anc1, dist = dist))
}

#' @name rsurv_prep
#' @export
rsurv_prep <- function(x){
  if(class(x) == "flexsurvreg"){
    d <- x$dlist
    param <- vector(length(d$pars) + 1, mode = "list")
    loc.indx <- which(d$pars == d$location)
    param[[1]] <- x$coef[c(x$basepars[loc.indx],
                           x$covpars[x$mx[[d$location]]])]
    if (length(d$pars) > 1){
      anc.indx <- which(d$pars != d$location)
      anc.names <- d$pars[anc.indx]
      for (i in 1:length(anc.names)){
        param[[i+1]] <- x$coef[c(x$basepars[anc.indx[i]],
                                 x$covpars[x$mx[[anc.names[i]]]])]
      }
    }
    param[[length(param)]] <- x$dlist$name
    if (x$dlist$name == "weibull.quiet"){
      param[[length(param)]] <- "weibull"
    }
    if (length(d$pars) <= 1){
      names.pars <- "location"
    } else{
      names.pars <-  c("location", paste0("anc", seq(1, length(param) - 2)))
    }
    names(param) <- c(names.pars, "dist")
  } else{
    print(paste0("rsurv_prep does not work for an object of class ", class(x)))
  }
  return(param)
}