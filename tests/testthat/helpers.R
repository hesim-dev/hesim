# functions from flexsurv package
bfp <- function (x, powers = c(1, 2)) {
  nobs <- length(x)
  npoly <- length(powers)
  X <- matrix(0, nrow = nobs, ncol = npoly)
  x1 <- ifelse(powers[1] != rep(0, nobs), x^powers[1], log(x))
  X[, 1] <- x1
  if (npoly >= 2) {
      for (i in 2:npoly) {
          if (powers[i] == powers[(i - 1)]) 
              x2 <- log(x) * x1
          else x2 <- ifelse(powers[i] != rep(0, nobs), x^powers[i], 
              log(x))
          X[, i] <- x2
          x1 <- x2
      }
  }
  X
}

