#' Simulate markov cohort model
#'
#' Simulate state transitions in a Markov cohort model.
#'
#' @param z0 Matrix of initial state vector (i.e. distribution of cohort
#' in each state during cycle 0) for each simulation.
#'
#' @param ncycles Number of cycles to run model.
#'
#' @param P Matrix with each row containing vector of transition probabilities for each
#' simulation and cycle. Vector of transition probabilities msut fill the transition matrix
#' rowwise. Also accepts a single vector of transition probabilities if transition matrix is
#'  constant accross simulations and cycles.
#'
#' @details The code is written in c++ to minmize simulation time.
#'
#' @return List with the following components.
#' \item{state}{Matrix containing number in cohort in each state during each cycle for each
#' simulation.}
#'
#' @export
markov_trans <- function(z0, ncycles, P, nsims = NULL){
  # z0 as a vector
  if (is.vector(z0)){
    z0 <- matrix(z0, nrow = 1)
  }
  if (nrow(z0) == 1) {
    if (is.null(nsims)){
      stop("nsims must be specified if z0 is a vector")
    }
    if (nsims > 1){
      z0 <- z0[rep(1, nsims), ]
    }
  }

  # z0 as a matrix
  if(!is.null(nsims) & nrow(z0) > 1){
    stop("Only specifiy nsims if z0 is constant accross simulations")
  }
  if (nrow(z0) > 1){
    nsims <- nrow(z0)
  }

  # P as a vector
  N <- nsims * ncycles
  if (is.vector(P)){
    P <- matrix(P, nrow = 1)
  }
  P <-  if (nrow(P) == 1 & N > 1) P[rep(1, N), ]

  # P as a matrix
  if (nrow(P) != N){
    stop("Number of rows in P times must equal number of simulations times number of cycles.")
  }

  # c++ function
  return(markov_transC(z0, ncycles, P))
}

#' @export
array_convert <- function(x){
  if (is.vector(x)){
    x <- array(x, dim = c(length(x), 1, 1))
  }
  if (is.matrix(x)){
    x <- array(x, dim = c(nrow(x), ncol(x), 1))
  }
  return(x)
}

