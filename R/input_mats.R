# Input matrices class ---------------------------------------------------------
#' Input matrices for a statistical model
#' 
#' @description 
#' An object of class `input_mats` contains input matrices
#' for simulating a statistical model. Consists of (i) input matrices, `X`, and 
#' (ii) [metadata][id_attributes()] used to index each matrix in `X`. 
#' 
#' Once created, an `input_mats` object can be converted 
#' to a [`data.table::data.table`] with `as.data.table()`, which is a helpful way to check that 
#' the object is as expected. The `print()` method summarizes the object and 
#' prints it to the console. 
#' 
#' More details are provided under "Details" below. 
#' 
#' 
#' @param X A list of input matrices for predicting the values of each parameter 
#' in a statistical model. May also be a list of lists of input matrices when a 
#' list of separate models is fit (e.g., with [flexsurvreg_list()]).
#' @param ... For `input_mats()`, arguments to pass to [id_attributes()]. For `print()`,
#' arguments to pass to [`data.table::print.data.table()`].
#' 
#' @details 
#' `input_mats` objects are used with [`params`] objects to simulate
#' disease models, cost models, and/or utility models. Each column of `$X` contains
#' variables from the `params` object and a given row corresponds to a combination
#' of the ID variables. An input matrix must always have rows for the treatment
#' strategies (`strategy_id`) and patients (`patient_id`); it may optionally 
#' have rows for health variables (`state_id` or `transition_id`) and time 
#' intervals (`time_id`). The rows must be sorted by prioritizing (i) `strategy_id`,
#' (ii) `patient_id`, (iii) the health related ID variable 
#' (either `state_id` or `transition_id`) and (iv) `time_id`.
#' 
#' While `input_mats` objects can be created directly with [input_mats()], it 
#' is rarely a good idea to do so. They are typically created as the 
#' `input_data` field when creating model objects (e.g., with 
#' [create_IndivCtstmTrans()], [create_CohortDtstmTrans()], and 
#' [create_PsmCurves()]). Internally, these functions 
#' create the input matrices using [create_input_mats()] methods, which ensure
#' that they are in the correct format. Users may also use [create_input_mats()] 
#' methods, but there is not usually a good reason to do so. 
#' 
#' `as.data.table.input_mats()` will convert input matrices into a single 
#' `data.table()` that column binds the ID variables and the unique combinations 
#' of variables contained in the elements of `$X`. `print.input_mats()` prints
#' a call to `as.data.table()` and provides additional information about the
#' ID variables.
#' 
#' @example man-roxygen/example-input_mats.R
#' @seealso See [IndivCtstmTrans()] and [PsmCurves()] for examples in which the
#' `input_data` field of an instance of a model class is an `input_mats` object.
#' 
#' @export
input_mats <- function(X, ...){
  object <- new_input_mats(X, ...)
  check(object)
  return(object)
}

new_input_mats <- function(X, ...){
  stopifnot(is.matrix(X) | is.list(X) | is.null(X)) 
  object <- c(list(X = X),
                   do.call("new_id_attributes", list(...)))
  object[sapply(object, is.null)] <- NULL
  class(object) <- "input_mats"
  return(object)
}

#' @rdname check
check.input_mats <- function(object, ...){
  # Check X
  if (!is.list(object$X)){
    stop("'X' must be a list or a list of lists.", call. = FALSE)
  }
  X <- flatten_lists(object$X)
  if(!all(sapply(X, function(y) inherits(y, "matrix")))){
    stop("'X' must be a list or list of lists of matrices.", call. = FALSE)
  }
  X_nrows <- sapply(X, nrow)
  if (is.list(object$X)){
    if (!all(X_nrows[1] == X_nrows)){
      stop("The number of rows in each matrix in 'X' must be the same.",
           call. = FALSE)
    }
  }
  
  # Check that rows in X are consistent with ID variables
  id_vars <- c("strategy_id", "patient_id", "state_id", "transition_id",
               "time_id")
  id_vars_n <- c("n_strategies", "n_patients", "n_states", "n_transitions",
                 "n_times")
  for (i in 1:length(id_vars)){
    if (!is.null(object[[id_vars[i]]])){
      if(length(object[[id_vars[i]]]) != X_nrows[1]){
        msg <- paste0("The length of '", id_vars[i], "' does not equal the number of rows in the ",
                      "'X' matrices.")
        stop(msg, call. = FALSE)
      }
    }
  }
  
  # Check ID attributes
  id_args <- object[names(object)[names(object) != "X"]]
  check(do.call("id_attributes", id_args))
}

# Convert input matrices to data tables ----------------------------------------
#' @param x An [`input_mats`] object.
#' @rdname input_mats
#' @export
as.data.table.input_mats <- function(x, ...) {
  
  # Get ID columns
  id_dt <- make_id_data_table(x)
  
  # Combine all X matrices
  xl <- lapply(flatten_lists(x$X), as.data.table)
  x_dt <- NULL
  for (i in 1:length(xl)) {
    cols_i <- colnames(xl[[i]])
    new_cols <- cols_i[!cols_i %in% colnames(x_dt)]
    if (length(new_cols) > 0) {
      x_dt <- cbind(x_dt, xl[[i]][, new_cols, with = FALSE])
    }
  }
  
  # Create a single data.table
  tbl <- cbind(id_dt, x_dt)
  setattr(tbl, "id_vars", attr(id_dt, "id_vars"))
  
  # Return
  tbl
}

# Print method for input matrices --------------------------------------------
#' @rdname input_mats
#' @export
print.input_mats <- function(x, ...) {
  x_dt <- as.data.table(x)
  id_vars <- attr(x_dt, "id_vars")
  size <- unlist(x[grepl("n_", names(x))])
  
  # Printing
  cat("An \"input_mats\" object \n\n")
  cat("Column binding the ID variables with all variables contained in the X matrices:\n")
  print(as.data.table(x), ...)
  cat("\n")
  cat("Number of unique values of ID variables:\n")
  print(size)
  cat("\n")
  if ("time_intervals" %in% names(x)) {
    cat("Time intervals:\n")
    print(x$time_intervals)
  }
  invisible(x)
}

# Helper functions to create input matrices ------------------------------------
size_id_map <- function(){
  c(strategy_id = "n_strategies", 
    patient_id = "n_patients",
    state_id = "n_states",
    transition_id = "n_transitions",
    time_id = "n_times")
}

get_input_mats_id_vars <- function(data){
  map <- size_id_map()
  res <- list() 
  id_vars <- attr(data, "id_vars")
  for (i in 1:length(id_vars)){
    res[[id_vars[i]]] <- data[[id_vars[i]]]
    res[[map[id_vars[i]]]] <- length(unique(data[[id_vars[i]]]))
    if (id_vars[[i]] == "time_id"){
      res[["time_intervals"]] <- attr(data, "time_intervals")
    }
  }
  if ("grp_id" %in% colnames(data)) res[["grp_id"]] <- data[["grp_id"]]
  if ("patient_wt" %in% colnames(data)) res[["patient_wt"]] <- data[["patient_wt"]]
  return(res)
}

#' Check input data argument for `create_input_mats`
#' 
#' Check that input data argument for `create_input_mats` exists and that it is
#' of the correct type. 
#' @param input_data An object of class "expanded_hesim_data" returned by the function
#'  [expand.hesim_data()]. 
#' @keywords internal
#' @return If all tests passed, returns nothing; otherwise, throws an exception.
check_input_data <- function(input_data){
  if (!inherits(input_data, "expanded_hesim_data")){
    stop("'input_data' must be of class 'expanded_hesim_data'.")
  } 
  if (!inherits(input_data, "data.table") & !inherits(input_data, "data.frame")){
    stop("'input_data' must inherit from either 'data.table' or 'data.frame'.")
  }   
  if (!inherits(input_data, "data.table")){
    setattr(input_data, "class", c("expanded_hesim_data", "data.table", "data.frame"))
  } 
}

extract_X <- function(coef_mat, data){
  varnames <- colnames(coef_mat)
  if (is.null(varnames)){
    stop("Variable names for coefficients cannot be NULL.",
         call. = FALSE)
  }
  if(!all(varnames %in% colnames(data))){
    stop("Not all variables in 'object' are contained in 'input_data'.",
         call. = FALSE)
  }  
  X <- as.matrix(data[, varnames, with = FALSE])
  if (!is.numeric(X)) {
    stop("'input_data' must only include numeric variables.",
         call. = FALSE)
  }
  return(X)
}

get_terms <- function(object){
  tt <- stats::terms(object)
  return(stats::delete.response(tt))
}

# Generic method for creating input matrices  ----------------------------------
#' Create input matrices
#' 
#' `create_input_mats()` is a generic function for creating an object of class
#' [`input_mats`]. Model matrices are constructed based on the 
#' variables specified in the model `object` and the data specified in `input_data`.
#' `create_input_mats()` is not typically called by users directly, but is 
#' instead used by functions that create model objects (e.g., 
#' [create_IndivCtstmTrans()], [create_CohortDtstmTrans()], 
#' [create_PsmCurves()]).
#' @param object An object of the appropriate class. 
#' @param input_data An object of class `expanded_hesim_data` returned by
#' [expand.hesim_data()]. It is used to look for the variables needed to create 
#' an input matrix for use in a statistical models and the ID variables for 
#' indexing rows in the input matrix. 
#' @param ... Further arguments passed to [model.matrix()].
#' @return An object of class `input_mats`.
#' @seealso [input_mats()]
#' @example man-roxygen/example-create_input_mats.R
#' @export
#' @keywords internal
#' @rdname create_input_mats
create_input_mats <- function (object, ...) {
  if (missing(object)){
    stop("'object' is missing with no default.")
  }
  UseMethod("create_input_mats", object)
}

# Create input matrices from formula  ------------------------------------------
formula_list_rec <- function(object, input_data, ...){
  x <- vector(mode = "list", length = length(object))
  names(x) <- names(object)
  for (i in 1:length(x)){
    if (inherits(object[[i]], "formula")){
      x[[i]] <- stats::model.matrix(object[[i]], data = input_data, ...)
    } else{
      x[[i]] <- formula_list_rec(object[[i]], data = input_data, ...)
    }
  }
  return(x)
}

#' Create input matrices from formula
#' 
#' This is an internal function for creating input matrices from formulas. It
#' is currently used in some unit tests.
#' @inheritParams create_input_mats
#' @inherit create_input_mats return
#' @export
#' @keywords internal
#' @seealso [create_input_mats()]
create_input_mats.formula_list <- function(object, input_data, ...){
  check_input_data(input_data)
  X_list <- formula_list_rec(object, input_data, ...)
  args <- c(list(X = X_list),
            get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

# Create input matrices from lm  -----------------------------------------------
#' @export 
#' @rdname create_input_mats
create_input_mats.lm <- function(object, input_data, ...){
  check_input_data(input_data)
  terms <- get_terms(object)
  X <- stats::model.matrix(terms, data = input_data, ...)
  args <- c(list(X = list(mu = X)),
           get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

# Create input matrices from flexsurvreg  --------------------------------------
create_input_mats_flexsurvreg_X <- function(object, input_data, ...){
  
  # Based on flexsurv:::form.model.matrix()
  mfo <- stats::model.frame(object)
  
  ## Error messages for missing variables in "input_data"
  covnames <- attr(mfo, "covnames")
  missing.covs <- unique(covnames[!covnames %in% names(input_data)])
  if (length(missing.covs) > 0){
    missing.covs <- sprintf("\"%s\"", missing.covs)
    plural <- if (length(missing.covs)>1) "s" else ""
    stop(sprintf("Value%s of covariate%s ",plural, plural),
         paste(missing.covs, collapse=", "), " not supplied in \"input_data\"")
  }
  
  ## as in predict.lm
  tt <- attr(mfo, "terms")
  Terms <- stats::delete.response(tt)
  mf <- stats::model.frame(Terms, input_data, xlev = stats::.getXlevels(tt, mfo))
  if (!is.null(cl <- attr(Terms, "dataClasses")))
    stats::.checkMFClasses(cl, mf)
  
  ## Return one model matrix for each parameter
  pars <- object$dlist$pars
  X_list <- vector(mode = "list", length = length(pars))
  names(X_list) <- pars
  for (i in 1:length(pars)){
    form <- object$all.formulae[[pars[i]]]
    if (is.null(form)){
      form <- stats::formula(~1)
    } else{
      form <- stats::delete.response(stats::terms(form))
    }
    X_list[[i]] <- stats::model.matrix(form, mf, ...)
  }
  return(X_list)
}

#' @export
#' @rdname create_input_mats
create_input_mats.flexsurvreg <- function(object, input_data,...){
  check_input_data(input_data)
  X_list <- create_input_mats_flexsurvreg_X(object, input_data, ...)
  args <- c(list(X = X_list),
           get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

#' @export
#' @rdname create_input_mats
create_input_mats.flexsurvreg_list <- function(object, input_data,...){
  check_input_data(input_data)
  X_list_2d <- vector(mode = "list", length = length(object))
  names(X_list_2d) <- names(object)
  for (i in 1:length(object)){
    X_list_2d[[i]] <- create_input_mats_flexsurvreg_X(object[[i]], input_data, ...)
  }
  args <- c(list(X = X_list_2d),
           get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

# Create input matrices from partsurvfit  --------------------------------------
#' @export
#' @rdname create_input_mats
create_input_mats.partsurvfit <- function(object, input_data, ...){
  check_input_data(input_data)
  return(create_input_mats.flexsurvreg_list(object$models, input_data, ...))
}

# Create input matrices from params_lm  ----------------------------------------
#' @export 
#' @rdname create_input_mats
create_input_mats.params_lm <- function(object, input_data, ...){
  check_input_data(input_data)
  X <- extract_X(object$coefs, input_data)
  args <- c(list(X = list(mu = X)),
            get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

# Create input matrices from params_surv  --------------------------------------
create_input_mats_params_surv_X <- function(object, input_data){
  X_list <- vector(mode = "list", length = length(object$coefs))
  names(X_list) <- names(object$coefs)
  for (i in 1:length(X_list)){
    X_list[[i]] <- extract_X(object$coefs[[i]], input_data)
  }
  return(X_list)
}

#' @export 
#' @rdname create_input_mats
create_input_mats.params_surv <- function(object, input_data, ...){
  check_input_data(input_data)
  X_list <- create_input_mats_params_surv_X(object, input_data)
  args <- c(list(X = X_list),
            get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

#' @export 
#' @rdname create_input_mats
create_input_mats.params_surv_list <- function(object, input_data, ...){
  X_list_2d <- vector(mode = "list", length = length(object))
  names(X_list_2d) <- names(object)
  for (i in 1:length(object)){
    X_list_2d[[i]] <- create_input_mats_params_surv_X(object[[i]], input_data)
  }
  args <- c(list(X = X_list_2d),
            get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

# Create input matrices from multinom  -----------------------------------------
create_input_mats_multinom_X <- function(object, input_data, ...){
  check_input_data(input_data)
  terms <- get_terms(object)
  if (!is.null(attr(terms(object), "offset"))){
    stop("An offset is not supported.", call. = FALSE)
  }
  m <- stats::model.frame(terms, input_data, na.action = stats::na.omit,
                           xlev = object$xlevels)
  if (!is.null(cl <- attr(terms, "dataClasses")))
    stats::.checkMFClasses(cl, m)
  return(stats::model.matrix(terms, m, contrasts = object$contrasts))
}

#' @export 
#' @rdname create_input_mats
create_input_mats.multinom <- function(object, input_data, ...){
  X <- create_input_mats_multinom_X(object, input_data, ...)
  args <- c(list(X = X),
            get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

#' @export 
#' @rdname create_input_mats
create_input_mats.multinom_list <- function(object, input_data, ...){
  X_list <- vector(mode = "list", length = length(object))
  names(X_list) <- names(object)
  for (i in 1:length(object)){
    X_list[[i]] <- create_input_mats_multinom_X(object[[i]], input_data, ...)
  }
  args <- c(list(X = X_list),
            get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

# Create input matrices from mlogit  -------------------------------------------
#' @export 
#' @rdname create_input_mats
create_input_mats.params_mlogit_list <- function(object, input_data, ...){
  X_list <- vector(mode = "list", length = length(object))
  for (i in 1:length(object)){
    X_list[[i]] <- extract_X(object[[i]]$coefs[, , 1], input_data)
  }
  args <- c(list(X = X_list),
            get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}
