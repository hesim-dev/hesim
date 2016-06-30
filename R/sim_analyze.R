#' Analyze simulated multi-state model
#'
#' Analyze a simulated multi-state model
#'
#'
#' @param sim A simulated multi-state model from \code{\link{sim_msm}}.
#'
#' @param r Discount rate for the model. Default is 0.03. Implies a discount factor of
#' \eqn{1/(1 + r)^t} where t is the cycle number.
#'
#' @param x Data matrix
#'
#' @param beta Matrix of regression coefficients corresponding to data matrix \code{x}.
#' One row for each state.
#'
#' @param poly.beta Matrix of coefficeints for piecewise time polynomials. One row for each
#' state.
#'
#' @param poly.deg Matrix indicating degree of polynomial on each interval. One row for each
#' state.
#'
#' @param knots Matrix indicating location of knots for piecewise polynomials. One row for
#' each state.
#' @return Vector of present value of quantity for each simulated
#' individual during time spent in state.
#'
#' @name sim_analyze
NULL

#' @rdname sim_analyze
#' @export
sim_pv <- function(sim, tmat, r = .03, x = NULL, beta, poly.beta = NULL, poly.deg = NULL,
                   knots = NULL, name = NULL, agevar = NULL){
  n.states <- length(unique(sim$state))
  zeros <- rep(0, n.states)
  if(is.data.frame(x)){
    x <- as.matrix(x)
  }
  if(is.null(knots)){
      knots <- cbind(zeros, zeros)
  } else{
      knots <- cbind(zeros, knots, zeros)
  }
  if (is.null(agevar)){
      agecol <- -1
      sim$age <- rep(0, nrow(sim))
  } else{
      agecol <- which(colnames(x) == agevar) - 1
  }
  if(!is.matrix(beta)){
    beta <- matrix(beta, nrow = n.states, ncol = 1)
  }
  if (is.null(x)){
    x <- matrix(1, nrow = length(unique(sim$id)), ncol = 1)
  }
  if(is.null(poly.beta)){
    poly.beta <- matrix(0, nrow = n.states, ncol = 1)
  }
  if(is.null(poly.deg)){
    poly.deg <- matrix(0, nrow = n.states, ncol = 1)
  }
  absorbing <- absorbing(tmat) - 1
  pv <- sim_msm_pvC(sim$id, sim$sim, sim$age, sim$state, sim$final, sim$time,
                    absorbing, r, x, agecol, beta, poly.beta, poly.deg, knots)
  colnames(pv) <- name
  return(pv)
}

#' @rdname sim_analyze
#' @export
sim_los <- function(sim, timevar){
  n <- length(unique(sim$id))
  los <- sim[, .(los = sum(get(timevar))), by = "state"]
  los$los <- los$los/n
  los <- los[order(state),]
  return(los)
}

#' @rdname sim_analyze
#' @export
sim_transprob <- function(sim, t){
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
