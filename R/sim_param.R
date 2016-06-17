#' Extract parameters for simulation from a multi-state model
#'
#' Extract parameters for simulation from a multi-state model fit using \code{flexsurvreg}
#' from the \code{flexsurv} package
#'
#' @param x A list of models fit with \code{flexsurvreg}
#'
#' @return List with the following parameters:
#' \item{loc_beta}{Regression coefficients for location parameter.}
#' \item{dist}{Distribution used for model fit.}
#' \item{par}{1st ancillary parameter in survvial model.}
#'
#' @export
sim_param <- function(x){
  loc_beta <- anc1 <- list()
  dist <- rep(NA, length(x))
  for (i in 1:length(x)){
    param <- sim_param1(x[[i]])
    loc_beta[[i]] <- param$location
    anc1[[i]] <- param$anc1
    dist[i] <- param$dist
  }
  loc_beta <- list_to_array(loc_beta)
  anc1 <- unlist(anc1)
  dist <- unlist(dist)
  return(list(loc_beta = loc_beta, anc1 = anc1, dist = dist))
}

sim_param1 <- function(x){
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
    print(paste0("SimPrep does not work for an object of class ", class(z)))
  }
  return(param)
}

