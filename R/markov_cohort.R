#' Simulate markov cohort model
#'
#' Simulate a (Bayesian) Markov cohort model.
#'
#' @param z0 Initial state vector (i.e. distribution of cohort in each state during cycle 0).
#'
#' @param ncycles Number of cycles to run model.
#'
#' @param discount Discount rate for the model. Default is 0.03. Implies a discount factor of
#' \eqn{1/(1 + discount)^t} where t is the cycle number.
#'
#' @param P An array of transition probabilities. Length of 1st dimension (# of rows) is
#' equal to the number of transition probability parameters. Transition probability parameters
#' should be ordered so that they fill a transition matrix by column. Length of 2nd dimension
#' (# of columns) is equal to the number of simulated (or user chosen) draws of the parameters.
#' Length of 3rd dimension is equal to the number of times the transition matrix changes over
#' time. Use \code{P_indx} to determine which array element should be used for each cycle.
#' Accepts a vector of transition probabilities if there is only 1 value for each parameter or a
#' matrix of transition probabilities if there are multiple simulated values of the parameters
#' but they do not vary over time.
#'
#' @param costs An arrary of costs for each state in model. Length of 1st dimension (rows) equal
#' to the number of model states; length of 2nd dimension (columns) equal to the number of
#' simulated (or user chosen) draws of the parameters; and length of the 3rd dimension
#' equal to the number of times costs change over time. Use \code{cost_indx} to determine which
#' array element should be used for each cycle. Accepts a vector of costs if there is only 1
#' value for each cost parameter or a matrix of costs if there are multiple simulated values of
#' the cost parameters but they do not vary over time.
#'
#' @param qol Arrary of quality of life weights for each state. Dimensions identical to
#' \code{costs} dimensions.
#'
#' @param P_indx Index vector of length ncycles. Each element corresponds to a cycle
#' and is equal to the element of the third dimension of the transition probability array,
#' \code{P}, that should be used for that cycle.
#'
#' @param cost_indx Index vector of length ncycles. Each element corresponds to a cycle
#' and is equal to the element of the third dimension of the cost array,
#' \code{costs}, that should be used for that cycle.
#'
#' @param qol_indx Index vector of length ncycles. Equivalent to \code{cost_indx}
#' but for quality of life.
#'
#' @details The code is written in c++ to minmize simulation time.
#'
#' @return List with the following components.
#' \item{state}{Matrix containing number in cohort in each state during each cycle.}
#' \item{costs}{Vector of total costs in each cycle.}
#' \item{qalys}{Vector of quality adjusted life years (QALYs) in each cycle}
#'
#' @export
markov_cohort <- function(z0, ncycles, discount = .03, nsims = 1,
                         P, costs, qol,
                         P_indx = NULL, cost_indx = NULL, qol_indx = NULL){
  # conditions for transition matrix
  P <- .arrayConvert(P)
  P_indx <-  if (dim(P)[3] == 1) rep(1, ncycles)

  # conditions for costs
  costs <- .arrayConvert(costs)
  cost_indx <- if (dim(costs)[3] == 1) rep(1, ncycles)

  # conditions for quality of life weights
  qol <- .arrayConvert(qol)
  qol_indx <- if (dim(qol)[3] == 1) rep(1, ncycles)

  # error messages
  if (!(ncycles) %in% c(length(P_indx), length(cost_indx), length(qol_indx))){
    stop("Index vector(s) must equal ncycles")
  }

  # c++ function
  return(markov_cohortC(z0, ncycles, discount, nsims,
                       P, costs, qol,
                       P_indx, cost_indx, qol_indx))
}

.arrayConvert <- function(x){
  if (is.vector(x)){
    x <- array(x, dim = c(length(x), 1, 1))
  }
  if (is.matrix(x)){
    x <- array(x, dim = c(nrow(x), ncol(x), 1))
  }
  return(x)
}

