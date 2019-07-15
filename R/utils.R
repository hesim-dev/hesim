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

check_dr <- function(dr){
  if(any(table(dr) > 1)){
    stop("You cannot specify the same discount rate twice.",
         call. = FALSE)
  }  
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

# R6 class for parameter tables (i.e., stateval_tbl, transprob_tbl)
ParamsTbl <- R6::R6Class("ParamsTbl",
  private = list(
    id_vars_all = c("sample", "strategy_id", "state_id", "transition_id",
                     "grp_id", "time_start"),
    id_vars = NULL,
    id_vars_msg = NULL,
    cols = NULL
  ),                       
                         
  public = list(
    tbl = NULL,
    dist = NULL,
    hesim_data = NULL, 
                
    initialize = function(tbl, dist, hesim_data){
      self$tbl <- tbl
      self$dist <- dist
      self$hesim_data <- hesim_data
      private$cols <- colnames(tbl)
    },
    
    add_time_intervals = function(tbl){
      if (!is.null(self$tbl$time_start)){
        time_intervals <- data.table(time_start = unique(self$tbl$time_start))
        time_intervals[, "time_stop" := shift(get("time_start"), type = "lead")]
        time_intervals[is.na(get("time_stop")), "time_stop" := Inf]
        time_intervals[, "time_id" := 1:nrow(time_intervals)]
        pos <- match(self$tbl$time_start, time_intervals$time_start)
        self$tbl[, "time_id" := time_intervals$time_id[pos]]
        self$tbl[, "time_stop" := time_intervals$time_stop[pos]]
        invisible(self)
      }
    },
    
    check_need_sample = function(){
      if (self$dist != "custom") {
        if ("sample" %in% colnames(self$tbl)){
          stop(paste0("If 'sample' is in 'tbl', then 'dist' must equal 'custom'."),
               call. = FALSE)
        }
      }      
    },
    
    check_dist = function(){
      if (self$dist == "norm"){
        if (!all(c("mean", "sd") %in% private$cols)){
          msg <- stop("If a normal distribution is specified, then tbl must ",
                      "contain the columns 'mean' and 'sd'.")
          stop(msg, call. = FALSE)         
        }
      } else if (self$dist == "beta"){
        if (!all(c("mean", "se") %in% private$cols) &
            !all(c("shape1", "shape2") %in% private$cols)){
          msg <- stop("If a beta distribution is specified, then tbl must either ",
                      "contain the columns 'mean' and 'se' or 'shape1' and 'shape2'.")
          stop(msg, call. = FALSE)      
        }
      } else if (self$dist == "gamma"){
        if (!all(c("mean", "se") %in% private$cols) &
            !all(c("shape", "rate") %in% private$cols) &
            !all(c("shape", "scale") %in% private$cols)){
          msg <- stop("If a gamma distribution is specified, then tbl must either ",
                      "contain the columns 'mean' and 'se', 'shape' and 'rate', ",
                      "or 'shape' and 'scale'.")
          stop(msg, call. = FALSE)        
        }
      } else if (self$dist == "lnorm"){
        if (!all(c("meanlog", "sdlog") %in% private$cols)){
          msg <- stop("If a lognormal distribution is specified, then tbl must ",
                      "contain the columns 'meanlog' and 'sdlog'.")
          stop(msg, call. = FALSE) 
        }
      } else if (self$dist == "unif"){
        if (!all(c("min", "max") %in% private$cols)){
          msg <- stop("If a uniform distribution is specified, then tbl must ",
                      "contain the columns 'min' and 'max'.")
          stop(msg, call. = FALSE)
        }
      } else if (self$dist == "fixed"){
        if (!all(c("est") %in% private$cols)){
          msg <- stop("If 'dist' = 'fixed', then tbl must ",
                      "contain the column 'est'.")
          stop(msg, call. = FALSE)
        }    
      }else if (self$dist == "dirichlet"){
        if (!all(c("alpha") %in% private$cols)){
          msg <- stop("If 'dist' = 'dirichlet', then tbl must ",
                      "contain the column 'alpha'.")
          stop(msg, call. = FALSE)
        }
      } else if (self$dist == "custom"){
        if (!all(c("sample", "value") %in% private$cols)){
          msg <- stop("If 'dist' = 'custom', then tbl must ",
                      "contain the columns 'sample' and 'value'.")
          stop(msg, call. = FALSE)
        }  
      }
    },
    
    check_need_hesim_data1 = function(var){
      if (is.null(self$tbl[[var]])){
        name <- switch(var,
                       "state_id" = "states",
                       "strategy_id" = "strategies",
                       "grp_id" = "patients")
        if (is.null(self$hesim_data[[name]])){
          msg <- paste0("If '", var, "' is not a column in 'tbl' ",
                        "then 'hesim_data' must be included as an argument ",
                        "and '",  name, "' must be an element of 'hesim_data'.")
          stop(msg, call. = FALSE)
        }
      }
      invisible(self)
    },
    
    check_need_hesim_data = function(vars){
      for (v in vars){
        self$check_need_hesim_data1(v)
      }
      invisible(self)
    },
    
    check_unique_rows = function(health_id){
      private$id_vars <- private$id_vars_all[which(private$id_vars_all %in% colnames(self$tbl))]
      if (length(private$id_vars) == 1){
        private$id_vars_msg <- private$id_vars
      } else{
        private$id_vars_msg <- paste0(paste(private$id_vars[1:length(private$id_vars) - 1], collapse = ", "),
                                      ", and ", private$id_vars[length(private$id_vars)])
      }
      if (!all(self$tbl[, .N, by = c(private$id_vars)]$N == 1)) {
        stop(paste0("There must only be one row for each combination of ",
                    private$id_vars_msg,
                    " in 'tbl'."),
             call. = FALSE)
      } 
      invisible(self)
    },
    
    get_unique_size = function(var){
      if (is.null(self$tbl[[var]])){
        n <- 1
      } else{
        n <- length(unique(self$tbl[[var]]))
      }
      return(n)  
    },
    
    check_size = function(){
      unique_sizes <- list()
      for (i in 1:length(private$id_vars_all)){
        unique_sizes[[i]] <- self$get_unique_size(private$id_vars_all[i])
      }
      expected_n <- prod(unlist(unique_sizes))
      if (nrow(self$tbl) != expected_n) {
        if (length(private$id_vars) == 1){
          stop(paste0("The number of rows in 'tbl' should equal ", expected_n, 
                      " which is the number of unique values of ",
                      private$id_vars_msg, " in 'tbl'")) 
        } else{
          stop(paste0("The number of rows in 'tbl' should equal ", expected_n, 
                      " which is the product of the number of unique values of ",
                      private$id_vars_msg, " in 'tbl'")) 
        }
      }
      invisible(self)
    },
    
    sort = function(){
      id_cols <- c("sample", "strategy_id", "patient_id", "grp_id", "state_id", "time_id",
                   "time_start", "time_stop") 
      pos <- which(id_cols %in% colnames(self$tbl))
      setcolorder(self$tbl, id_cols[pos])
      invisible(self)
    },
    
    set_attributes = function(health_id){
      if (health_id == "state_id"){
        setattr(self$tbl, "class", c("stateval_tbl", "data.table", "data.frame"))
        setattr(self$tbl, "state_id", self$hesim_data$states$state_id)
        setattr(self$tbl, "strategy_id", self$hesim_data$strategies$strategy_id)
      } else{
        setattr(self$tbl, "class", c("transprob_tbl", "data.table", "data.frame"))
      }
      setattr(self$tbl, "dist", self$dist)
      setattr(self$tbl, "patients", data.table(self$hesim_data$patients))
    },
    
    create_tbl = function(tbl = c("stateval_tbl", "transprob_tbl")){
      if (tbl == "stateval_tbl"){
        need_hesim_data_vars <- c("state_id", "strategy_id", "grp_id")
        health_id <- "state_id"
      } else{
        need_hesim_data_vars <- "grp_id"
        health_id <- "transition_id"
      }
      
      # Check
      self$check_dist()
      self$check_need_sample()
      self$check_need_hesim_data(need_hesim_data_vars)
      self$check_unique_rows()
      self$check_size()    
      
      # Modify table
      self$add_time_intervals()
      self$sort()
      self$set_attributes(health_id)      
      
    }
  )
)
