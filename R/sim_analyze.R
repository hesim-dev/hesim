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
sim_pv <- function(sim, r = .03, x, beta, poly.beta = NULL, poly.deg = NULL,
                   knots = NULL){
  n.states <- length(unique(sim.cea$state))
  zeros <- rep(0, n.states)
  if(is.null(knots)){
      knots <- cbind(zeros, zeros)
  } else{
      knots <- cbind(zeros, knots, zeros)
  }
  if ("age" %in% colnames(sim)){
      agecol <- which(colnames(sim) == "age") - 1
  } else{
    agecol <- -1
    sim$age <- rep(0, nrow(sim))
  }
  if(!is.matrix(beta)){
    beta <- matrix(beta, nrow = n.states, ncol = 1)
  }
  if(is.null(poly.beta)){
    poly.beta <- matrix(0, nrow = n.states, ncol = 1)
  }
  if(is.null(poly.deg)){
    poly.deg <- matrix(0, nrow = n.states, ncol = 1)
  }
  pv <- sim_msm_pvC(sim$id, sim$sim, sim$age, sim$state, sim$time,
                    r, x, agecol, beta, poly.beta, poly.deg, knots)
  return(pv)
}

#' @rdname sim_analyze
#' @export
sim_los <- function(sim){
  s = copy(sim)
  s[, n := seq(1, .N), by = "id"]
  s[, N := .N, by = "id"]
  s.last <- s[n == N]
  return(mean(s.last$time))
}
