#' Length of stay 
#'
#' Length of stay from continuous time simulation.
#'
#'
#' @param sim An object of class \code{ctsim}.
#' @return Length of stay for each sampled parameter set. 
#'
#' @export
los_ctsim <- function(sim, timevar){
  n <- length(unique(sim$id))
  los <- sim[, .(los = sum(get(timevar))), by = "state"]
  los$los <- los$los/n
  los <- los[order(state),]
  return(los)
}

#' Transition probabilities 
#'
#' Transition probabilities for continuous time simulation.
#'
#'
#' @param sim An object of class \code{ctsim}.
#' @return Length of stay for each sampled parameter set. 
#'
#' @export
transprob_ctsim <- function(sim, t){
  n.states <- length(unique(sim$state))
  simindivs <- length(unique(sim$id)) * length(unique(sim$sim))
  transprob <- sim_transprobC(sim$state, sim$time, sim$final, t, simindivs, n.states)
  transprob <- matrix(transprob, nrow = length(t), ncol = n.states, byrow = T)
  if (t[1] == 0){
    start.state <- sim$state[1] + 1
    transprob[1, start.state] <- 1
    transprob[1, -start.state] <- 0
  }
  return(transprob)
}
