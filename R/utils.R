# Check ------------------------------------------------------------------------
#' Input validation for class objects
#' 
#' \code{check} is a generic function for validating the inputs of class objects.
#' @param object object to check.
#' @param inner_class When checking a list of objects, the class of elements within
#' the inner most list.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return If validation is successful, returns the object in question; otherwise,
#' informs the user that an error has occurred.  
#' @keywords internal
check <- function (object, ...) {
  UseMethod("check")
}

# Additional utility methods ---------------------------------------------------
#' Absorbing states
#' 
#' Returns a vector of absorbing states from a transition matrix.
#' @param trans_mat A transition matrix in the format from the \link[mstate]{mstate} package. 
#' See \link{IndivCtstmTrans}.
#' @keywords internal
absorbing <- function(trans_mat){
  which(apply(trans_mat, 1, function(x) all(is.na(x))))
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

# Get the object containing ID attributes
get_id_object <- function(x){
  if (is.null(x$input_mats)){
    return(x$params)
  } else{
    return(x$input_mats)
  }
}

# R6 class for parameter tables (i.e., stateval_tbl, transprob_tbl) ------------
ParamsTbl <- R6::R6Class("ParamsTbl",
  private = list(
    id_vars_all = NULL,
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
    
    set_id_vars_all = function(){
      if ("patient_id" %in% colnames(self$tbl)){
        private$id_vars_all <-  c("sample", "strategy_id", "state_id", "transition_id",
                                  "patient_id", "time_start")
      } else {
        private$id_vars_all <-  c("sample", "strategy_id", "state_id", "transition_id",
                                  "grp_id", "time_start")
      } 
    },
    
    add_time_intervals = function(){
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
      } else if (self$dist == "dirichlet"){
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
      if (var %in% c("patient_id", "grp_id")){ # Patient ID or Group ID
        if (is.null(self$tbl[["patient_id"]]) & is.null(self$tbl[["grp_id"]])){
          if (is.null(self$hesim_data[["patients"]])){
            msg <- paste0("If 'either 'patient_id' or 'grp_id' is not a column in 'tbl' ",
                          "then 'hesim_data' must be included as an argument ",
                          "and 'patients' must be an element of 'hesim_data'.")
            stop(msg, call. = FALSE)
          }          
        } 
      } else { # Other (strategy ID and state ID)
        if (!var %in% c("grp_id", "patient_id") & is.null(self$tbl[[var]])) {
          name <- switch(var,
                         "state_id" = "states",
                         "strategy_id" = "strategies")
          if (is.null(self$hesim_data[[name]])){
            msg <- paste0("If '", var, "' is not a column in 'tbl' ",
                          "then 'hesim_data' must be included as an argument ",
                          "and '",  name, "' must be an element of 'hesim_data'.")
            stop(msg, call. = FALSE)
          }
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
      self$set_id_vars_all()
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

# R6 class for creating economic model component from a parameter table --------
CreateFromParamsTbl <- R6::R6Class("CreateFromParamsTbl",
  public = list(
    object = NULL,
    n = NULL, 
    values = NULL,
    id_tbl = NULL,
    params = NULL,
    n_trans = NULL, # For transprob_tbl only
    n_states = NULL, # For transprob_tbl only
                    
    initialize = function(object, n){
      self$object <- object
      self$n <- n
      if (inherits(self$object, "transprob_tbl")){
        self$n_trans <- length(unique(self$object$transition_id))
        self$n_states <- sqrt(self$n_trans)
      }
    },
    
    random = function(){
      n_rows <- nrow(self$object)
      if (attr(self$object, "dist") == "norm"){
        self$values <- stats::rnorm(self$n * n_rows, 
                                    mean = self$object$mean, sd = self$object$sd)
      } else if (attr(self$object, "dist") == "beta"){
        if (all(c("shape1", "shape2") %in% colnames(self$object))){
          self$values <- stats::rbeta(self$n * n_rows, 
                                      shape1 = self$object$shape1, 
                                      shape2 = self$object$shape2)
        } else if (all(c("mean", "se") %in% colnames(self$object))){
          mom_params <- mom_beta(self$object$mean, self$object$se)
          self$values <- stats::rbeta(self$n * n_rows, 
                                      shape1 = mom_params$shape1, 
                                      shape2 = mom_params$shape2) 
        } 
      } else if (attr(self$object, "dist") == "gamma"){
        if (all(c("shape", "rate") %in% colnames(self$object))){
          self$values <- stats::rgamma(self$n * n_rows, 
                                       shape = self$object$shape, 
                                       rate = self$object$rate)
        } else if (all(c("shape", "scale") %in% colnames(self$object))){
          self$values <- stats::rgamma(self$n * n_rows, 
                                       shape = self$object$shape, 
                                       scale = self$object$scale)
        } else if (all(c("mean", "se") %in% colnames(self$object))){
          mom_params <- mom_gamma(self$object$mean, self$object$se)
          self$values <- stats::rgamma(self$n * n_rows, 
                                       shape = mom_params$shape, 
                                       scale = mom_params$scale) 
        } 
      } else if (attr(self$object, "dist") == "lnorm"){
        self$values <- stats::rlnorm(self$n * n_rows, 
                                     meanlog = self$object$meanlog, 
                                     sdlog = selfdevt$object$sdlog)
      } else if (attr(self$object, "dist") == "unif"){
        self$values <- stats::runif(self$n * n_rows, 
                                     min = self$object$min, 
                                     max = self$object$max) 
      } else if (attr(self$object, "dist") == "dirichlet"){
        alpha <- matrix(self$object$alpha, ncol = self$n_states, byrow = TRUE)
        self$values <- rdirichlet_mat(self$n, alpha)
      } else if (attr(self$object, "dist") == "fixed"){
        self$values <- rep(self$object$est, times = self$n)
      } else if (attr(self$object, "dist") == "custom"){
        self$values <- self$object$value
      }      
      invisible(self)
    },
    
    transform = function(){
      if (attr(self$object, "dist") == "dirichlet"){ # Dirichlet, only for transprob_tbl
        self$values <- aperm(array(c(aperm(self$values, perm = c(2, 1, 3))),
                                   dim = c(self$n_states, self$n_states, 
                                           dim(self$values)[3] * dim(self$values)[1]/self$n_states)),
                             perm = c(2, 1, 3))
      } else if (attr(self$object, "dist") == "custom"){ # Custom distribution
        setorderv(self$object, "sample") 
        n_samples <- length(unique(self$object$sample))
        if(inherits(self$object, "stateval_tbl")){
          self$values <- matrix(self$values, ncol = n_samples, byrow = FALSE)
        } else{
          self$values <- aperm(array(self$values,
                                     dim = c(self$n_states, self$n_states, 
                                             nrow(self$object)/self$n_trans)),
                               perm = c(2, 1, 3))
        }
        if (self$n < n_samples){
          samples <- sample.int(n_samples, self$n, replace = FALSE) 
        } else if (self$n > n_samples) {
          warning("'n' is larger than the number of unique values of 'sample' in 'object'.")
          samples <- sample.int(n_samples, self$n, replace = TRUE) 
        }
        if (self$n != n_samples){
          self$values <- self$values[, samples, drop = FALSE]
        }        
      } else{ # All other distributions
        if(inherits(self$object, "stateval_tbl")){
          self$values <- matrix(self$values, ncol = self$n, byrow = FALSE) 
        } else{
          self$values <- aperm(array(self$values,
                                     dim = c(self$n_states, self$n_states, 
                                             length(self$values)/self$n_trans)),
                               perm = c(2, 1, 3))
        }        
      }
      invisible(self)
    },
    
    expand = function(){
      if (attr(self$object, "dist") == "custom"){
        self$id_tbl <- self$object[sample == 1]
      } else{
        self$id_tbl <- copy(self$object)
      }
      if(inherits(self$object, "transprob_tbl")){
        self$id_tbl <- self$id_tbl[transition_id == 1]
      }          
      self$id_tbl[, ("obs_num") := 1:.N]  
      
      # Expand by strategy_id and/or state_id and time interval
      id_tbl_list <- list()
      id_vars <- c("strategy_id", "state_id", "time_start")
      i <- 1
      for (var in id_vars){
        if (is.null(self$id_tbl[[var]])){
          if (!is.null(attr(self$object, var))){
            id_tbl_i <- data.frame(tmp_var = attr(self$object, var))
            setnames(id_tbl_i, "tmp_var", var)
            id_tbl_list[[i]] <- id_tbl_i
            i <- i + 1
          }
        }
      }
      id_tbl_list <- c(list(data.frame(self$id_tbl)), id_tbl_list)
      self$id_tbl <- Reduce(function(...) merge(..., by = NULL), id_tbl_list)
      self$id_tbl <- data.table(self$id_tbl)
      
      # Expand by patient 
      if (is.null(self$id_tbl$patient_id)){
        merge <- TRUE
        if (is.null(self$id_tbl$grp_id)){ # If group ID is not specified
          self$id_tbl[, ("grp_id") := 1]
          patient_lookup <- data.table(patient_id = attr(self$object, "patients")$patient_id, 
                                       grp_id = 1)
        } else { # Else if group ID is specified
          if (is.null(attr(self$object, "patients")$grp_id)) { # If the patient lookup table does not exist
            setnames(self$id_tbl, "grp_id", "patient_id")
            merge <- FALSE
          } else{
            patient_lookup <- attr(self$object, "patients")[, c("patient_id", "grp_id"), 
                                                            with = FALSE] 
          }
        }
        if (merge){
          self$id_tbl <- merge(self$id_tbl, patient_lookup, by = c("grp_id"), allow.cartesian = TRUE,
                               sort = FALSE) 
        }  
      }
      
      # Sort
      health_id <- switch(class(self$object)[1],
                          stateval_tbl = "state_id",
                          transprob_tbl = "transition_id")
      if (is.null(self$id_tbl[["time_start"]])){
        setorderv(self$id_tbl, cols = c("strategy_id", "patient_id", health_id)) 
      } else{
        setorderv(self$id_tbl, cols = c("strategy_id", "patient_id", health_id, 
                                   "time_id")) 
      }
      
      # Expanded values
      if (inherits(self$object, "stateval_tbl")){
        self$values <- self$values[self$id_tbl$obs_num, , drop = FALSE]
      } else{
        self$values <- self$values[,, rep(self$id_tbl$obs_num, each = self$n)]
      }
      invisible(self)
    },
    
    create_params = function(){
      if (!is.null(self$object$time_id)){
        time_intervals <- unique(self$object[, c("time_id", "time_start", "time_stop")])
      } else{
        time_intervals <- NULL
      }
      if (is.null(self$id_tbl$state_id)){
        n_states <- NULL
      } else{
        n_states <-  length(unique(self$id_tbl$state_id))
      }
      n_strategies <- length(unique(self$id_tbl$strategy_id))
      n_patients <- length(unique(self$id_tbl$patient_id))
      n_times <- nrow(time_intervals)
      
      if (inherits(self$object, "stateval_tbl")){
        self$params <- new_tparams_mean(value = self$values,
                                        n_samples = self$n,
                                        strategy_id = self$id_tbl$strategy_id,
                                        n_strategies = n_strategies,
                                        patient_id = self$id_tbl$patient_id,
                                        n_patients = n_patients,
                                        state_id = self$id_tbl$state_id,
                                        n_states = n_states,
                                        time_id = self$id_tbl$time_id,
                                        time_intervals = time_intervals,
                                        n_times = nrow(time_intervals))
      } else{
        self$params <- new_tparams_transprobs(value = self$values, 
                                              sample = rep(1:self$n,
                                                           each = dim(self$values)[3]/self$n),
                                              n_samples = self$n,
                                              strategy_id = rep(self$id_tbl$strategy_id, self$n),
                                              n_strategies = n_strategies,
                                              patient_id = rep(self$id_tbl$patient_id, self$n),
                                              n_patients = n_patients,
                                              state_id = rep(self$id_tbl$state_id, self$n),
                                              n_states = n_states,
                                              time_id = rep(self$id_tbl$time_id, self$n),
                                              time_intervals = time_intervals,
                                              n_times = n_times)
      }
    },
    
    prep = function(...){
      self$random()
      self$transform()
      self$expand()
      self$create_params()
    }
  )
)