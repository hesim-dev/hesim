#' Random generation for categorical distribution
#'
#' Draw random samples from a categorical distribution given a matrix of probabilities.
#'  \code{rcat} is vectorized and written in C++ for speed.  
#'
#' @param n Number of random observations to draw.
#' @param prob A matrix of probabilities where rows correspond to observations 
#' and columns correspond to categories.
#' @return A vector of random samples from the categorical distribution. The length of the sample is 
#' determined by n. The numerical arguments other than n are recycled so that the number of samples is 
#' equal to n.
#' @name rcat
#' @examples
#' p <- c(.2, .5, .3)
#' n <- 10000
#' pmat <- matrix(rep(p, n), nrow = n, ncol = length(p), byrow = TRUE)
#'
#' # rcat
#' set.seed(100)
#' ptm <- proc.time()
#' samp1 <- rcat(n, pmat)
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
rcat <- function(n, prob){
  if(!is.matrix(prob) & !is.vector(prob)){
      stop("prob must be a vector or a matrix.")
 } 
  if (is.vector(prob)){
    prob <- matrix(prob, nrow = 1)
  }
  if (n <= 0){
    stop("n must be greater than 0")
  }
  sim <- C_rcat_vec(n, prob) + 1
  return(as.numeric(sim))
}

#' Random generation for piecewise exponential distribution
#'
#' Draw random samples from an exponential distribution with piecewise rates.
#'  \code{rpwexp} is vectorized and written in C++ for speed.  
#'
#' @param n Number of random observations to draw.
#' @param rate A matrix of rates where rows correspond to observations 
#' and columns correspond to rates during specified time intervals. 
#' @param time A vector equal to the number of columns in \code{rate} giving the
#' times at which the rate changes
#' @return A vector of random samples from the piecewise exponential distribution. The length of the sample is 
#' determined by n. The numerical arguments other than n are recycled so that the number of samples is 
#' equal to n.
#' @name rpwexp
#' @examples
#' rate <- c(.6, 1.2, 1.3)
#' n <- 100000
#' ratemat <- matrix(rep(rate, n/2), nrow = n, 
#'                   ncol = 3, byrow = TRUE)
#' t <- c(0, 10, 15) 
#' ptm <- proc.time()
#' samp <- rpwexp(n, ratemat, t)
#' proc.time() - ptm
#' summary(samp)
#' @export
rpwexp <- function(n, rate = 1, time = 0){
  if(!is.matrix(rate) & !is.vector(rate)){
    stop("rate must be a vector or a matrix.")
  } 
  if (!is.vector(time)){
    stop("time must be a vector")
  }
  if (is.vector(rate)){
    if (length(rate) != length(time)){
      stop("length of time must be equal to the lenght of rate.")
    }
    rate <- matrix(rate, nrow = 1)
  }
  if (n <= 0){
    stop("n must be greater than 0")
  }
  if (ncol(rate) != length(time)){
    stop("length of time must be equal to the number of columns in rate.")
  }
  return(C_rpwexp_vec(n, rate, time))
}

#' Random generation for multiple Dirichlet distributions
#'
#' Draw random samples from multiple Dirichlet distributions.
#'\code{rdirichlet_mat} is vectorized and written in C++ for speed.  
#'
#' @param n Number of samples to draw.
#' @param alpha A matrix where each row is a separate vector of shape parameters.
#' @name rdirichlet_mat
#' @examples
#' alpha <- matrix(c(100, 200, 500, 50, 70, 75), ncol = 3, nrow = 2, byrow = TRUE)
#' samp <- rdirichlet_mat(100, alpha)
#' print(samp[, , 1:2])
#' @details This function is particularly useful for representing the distribution of 
#' transition probabilities in a transition matrix.
#' @return An array of matrices where each row of each matrix is a sample from the Dirichlet distribution.
#' @export
rdirichlet_mat <- function(n, alpha){
  if (n <= 0){
    stop("n must be greater than 0")
  }
  if (!is.matrix(alpha) & !is.vector(alpha)){
    stop("alpha must be a vector or a matrix")
  }
  if (is.vector(alpha)){
    alpha <- matrix(alpha, nrow = 1)
  }
  samp <- C_rdirichlet_mat(n, alpha)
  return(samp)
}

#' Random generation for generalized gamma distribution
#'
#' Draw random samples from a generalized gamma distribution using the 
#' parameterization from \code{flexsurv}. Written in C++
#' for speed. Equivalent to \code{flexsurv::rgengamma}.
#'
#' @param n Number of random observations to draw.
#' @param mu Vector of location parameters. 
#' and columns correspond to rates during specified time intervals. 
#' @param sigma Vector of scale parameters as described in \code{flexsurv}.
#' @param Q Vector of shape parameters.
#' @name fast_rgengamma
#' @examples
#' n <- 1000
#' m <- 2 ; s <- 1.7; q <- 1
#' ptm <- proc.time()
#' r1 <- fast_rgengamma(n, mu = m, sigma = s, Q = q)
#' proc.time() - ptm
#' ptm <- proc.time()
#' library("flexsurv")
#' r2 <- flexsurv::rgengamma(n, mu = m, sigma = s, Q = q)
#' proc.time() - ptm
#' summary(r1)
#' summary(r2)
#' 
#' @return A vector of random samples from the generalized gamma distribution. The length of the sample is 
#' determined by n. The numerical arguments other than n are recycled so that the number of samples is 
#' equal to n.
#' @export
fast_rgengamma <- function(n, mu = 0, sigma = 1, Q){
  return(C_rgengamma_vec(n, mu, sigma, Q))
}
