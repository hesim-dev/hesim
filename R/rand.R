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
  sim <- C_rcat(n, prob) + 1
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
  return(C_rpwexp(n, rate, time))
}

#' Random generation for multiple Dirichlet distributions
#'
#' Draw random samples from multiple Dirichlet distributions for use in 
#' transition probability matrices.
#' 
#'
#' @param n Number of samples to draw.
#' @param alpha A matrix where each row is a separate vector of shape parameters.
#' @param output Whether output should be an array or a matrix. 
#' @name rdirichlet_mat
#' @examples
#' alpha <- matrix(c(100, 200, 500, 50, 70, 75), ncol = 3, nrow = 2, byrow = TRUE)
#' samp <- rdirichlet_mat(100, alpha)
#' print(samp[, , 1:2])
#' @details This function is meant for representing the distribution of 
#' transition probabilities in a transition matrix. The `(i,j)` element of
#' `alpha` is a transition from state `i` to state `j`. The number of rows in
#' `alpha` must be less than or equal to the number of columns.
#'  \code{rdirichlet_mat} is vectorized and written in C++ for speed. 
#' @return If `output = "array"`, then an array of matrices is returned 
#' where each row of each matrix is a sample from the Dirichlet distribution.
#' If `output = "matrix"`, then a matrix is returned where each row contains
#' all elements of the sampled matrix from the Dirichlet distribution ordered rowwise; 
#' when `output = "data.frame"` or `output = "data.table"` the returned
#' value is a `data.frame` or `data.table` in the same format as the
#' matrix returned with `output = "matrix"`. 
#' @export
rdirichlet_mat <- function(n, alpha, output = c("array", "matrix", "data.frame", 
                                                "data.table")){
  output <- match.arg(output)
  if (n <= 0){
    stop("n must be greater than 0")
  }
  if (!(is.numeric(alpha) & length(dim(alpha)) <= 2)){
    stop("alpha must be a numeric matrix or vector")
  }
  if (!is.matrix(alpha)){
    alpha <- matrix(alpha, nrow = 1)
  }
  if (nrow(alpha) > ncol(alpha)){
    stop(paste0("The number of rows of 'alpha' must be less than or equal to ",
                "the number of columns."))
  }
  samp <- C_rdirichlet_mat(n, alpha)
  if (output %in% c("matrix", "data.frame", "data.table")){
    samp <- matrix(c(aperm(samp, perm = c(2, 1, 3))),
                   ncol = dim(samp)[1] * dim(samp)[2], byrow = TRUE)
    colnames(samp) <- paste0("prob_", 1:ncol(samp))
  }
  if (output == "data.frame"){
    samp <- data.frame(samp)
  } 
  if (output == "data.table"){
    samp <- data.table(samp)
  }
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
  return(C_rgengamma(n, mu, sigma, Q))
}
