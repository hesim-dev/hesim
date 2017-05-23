#' Simulate continuous time semi-markov (clock-reset) multi-state model
#'
#' Simulate continuous time semi-markov (clock-reset) multi-state model
#'
#' @param loc_beta Regression coefficients for location parameter. Array with 1st dimension (rows)
#' indexing random draw of coefficients, 2nd dimension (columns) indexing coefficients, and
#' 3rd dimension (slice of cube) equal to number of unique transitions.
#'
#' @param x Data matrix for coefficients. Number of columns must equal length of
#' 2nd dimension in loc_beta.
#'
#' @param dist Character vector indicating the probability distributions used for each
#'  transition. These include "weibull", "exponential", and "gompertz".
#'
#' @param tmat Matrix indicating model transitions.
#'
#' @param anc1 1st ancillary parameter (those other than the location parameter) in the model.
#'
#' @param maxt Time to simulate model until.
#'
#' @param maxage Maximum age in simulation.
#'
#' @param agevar Name of age variable in x.
#'
#' @details The code is written in c++ to minmize simulation time.
#'
#' @return Dataframe with the following columns:
#' \item{id}{Identification number for each individual.}
#' \item{sim}{Simulation number (for each draw of parameters).}
#' \item{state}{Model state during each transition.}
#' \item{time}{Time state was reached.}
#'
#' @export
#' @keywords internal
sim_ctsm <- function(loc_beta, x, dist, tmat, anc1, maxt, agevar = NULL, maxage = 101){
  loc_beta <- list_to_array(loc_beta)
  if (is.data.frame(x)){
    x <- as.matrix(x)
  }
  if (is.null(anc1)){
    anc1 = -1
  } else{
    anc1 = list_to_array(anc1)
  }
  if (is.null(agevar)){
      agecol <- -1
  } else{
      agecol <- which(colnames(x) == agevar) - 1
  }
  if (sum(!dist %in% c("exponential", "weibull", "gompertz")) > 0){
    stop("Distribution not recognized")
  }
  absorbing <- absorbing(tmat) - 1
  tmat[is.na(tmat)] <- 0
  sim <- sim_ctsmC(loc_beta, x, dist, tmat, anc1, absorbing, maxt, maxage, agecol)
  sim <- as.data.frame(sim)
  if (is.null(agevar)){
      colnames(sim) <- c("id", "sim", "state", "final", "time")
  } else{
    colnames(sim) <- c("id", "sim", "age", "state", "final", "time")
  }
  return(data.table(sim))
}


