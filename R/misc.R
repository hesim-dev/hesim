list_to_array <- function(L){
  if (is.matrix(L[[1]]) == T){
      array(unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)))
  } else if (is.vector(L[[1]]) == T){
      array(unlist(L), dim = c(1, length(L[[1]]), length(L)))
  } else{
      stop("List must contain matrices or vectors")
  }
}
