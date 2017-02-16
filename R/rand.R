#' Fast random number generation for categorical distribution
#'
#' Generate random numbers from a categorical distribution given a matrix of probabilities.
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