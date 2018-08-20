#' Absorbing states
#' 
#' Returns a vector of absorbing states from a transition matrix.
#' @param trans_mat A transition matrix in the format from the \link[mstate]{mstate} package. 
#' See \link{IndivCtstmTrans}.
#' @keywords internal
absorbing <- function(trans_mat){
  which(apply(trans_mat, 1, function(x) all(is.na(x))))
}

#' Input validation for class objects
#' 
#' \code{check} is a generic function for validating the inputs of class objects.
#' @param object object to check.
#' @param inner_class When checking a list of objects, the class of elements within
#' the inner most list.
#' @param ... Further arguments passed to or from other methods. Currently unused.
#' 
#' @return If validation is successful, returns the object in question; otherwise,
#' informs the user that an error has occurred.  
#' @keywords internal
check <- function (object, ...) {
  UseMethod("check")
}

#' Form a list from \code{...}
#' 
#' Form a list of objects from \code{...}.
#' @param ... Objects used to form a list.
#' @return A list of objects from \code{...}.
#' @keywords internal
create_object_list <- function(...){
  objects <- list(...)
  if(length(objects) == 1 & inherits(objects[[1]], "list")){
    objects <- objects[[1]]
  }
  return(objects)
}

# Create list of objects
check_object_list <- function(x, inner_class){
  for (i in 1:length(x)){
    if(!inherits(x[[i]], inner_class)){
      msg <- paste0("Each element in ... must be of class '", inner_class, "'")
      stop(msg, call. = FALSE)
    }
  } 
  return(x)
}

new_object_list <- function(..., new_class){
  objects <- create_object_list(...)
  class(objects) <- new_class
  return(objects)
}

object_list <- function(..., inner_class, new_class){
  res <- new_object_list(..., new_class = new_class)
  check_object_list(res, inner_class)
}

# Join objects at specified time points
check_joined_object <- function(x, inner_class, model_list){
  check_object_list(x$models, inner_class)
  
  if(model_list == FALSE){
     check_joined_times(x$models, x$times)
  } else {
    if(!is.list(x$times)){
      stop("'times' must be a list.", call. = FALSE)
    }
    for (i in 1:length(x$times)){
      check_joined_times(x$models[[i]], x$times[[i]])
    }
  } 
  return(x)
}

new_joined_object <- function(..., times, new_class){
  objects <- create_object_list(...)
  res <- list(models = objects, times = times)
  class(res) <- new_class
  return(res)
}


joined_object <- function(..., times, inner_class, new_class, model_list = FALSE){
  res <- new_joined_object(..., times = times, new_class = new_class)
  check_joined_object(res, inner_class, model_list)
}

check_joined_times <- function(objects, times){
  stopifnot(is.vector(times))
  stopifnot(is.numeric(times))
  stopifnot(!is.unsorted(times))
  if(length(objects) != (length(times) + 1)){
    stop("Length of joined models must equal 'times' + 1.",
         call. = FALSE)
  }
}


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

# List depth
list_depth <- function(list) {
  ifelse(is.list(list), 1L + max(sapply(list, list_depth)), 0L)
}

# Flatten a nested list
flatten_lists <- function(x) {
  if (!inherits(x, "list")) return(list(x))
  else return(unlist(c(lapply(x, flatten_lists)), recursive = FALSE))
}

