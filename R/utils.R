# list to array
list_to_array <- function(L){
  if (is.matrix(L[[1]]) == TRUE){
      array(unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)))
  } else if (is.vector(L[[1]]) == TRUE){
      array(unlist(L), dim = c(1, length(L[[1]]), length(L)))
  } else{
      stop("List must contain matrices or vectors")
  }
}

# Return vector of absorbing states
absorbing <- function(tmat){
  which(apply(tmat, 1, function(x) all(is.na(x))))
}
