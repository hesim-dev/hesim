# input_mats class -------------------------------------------------------------
#' Input matrices for a statistical model
#' 
#' Create an object of class `input_mats`, which contains inputs matrices
#' for simulating a statistical model. Consists of (i) input matrices, `X`, and 
#' (ii) [metadata][id_attributes()] used to index each matrix in `X`. 
#' More details are provided under "Details" below. 
#' 
#' @param X A list of input matrices for predicting the values of each parameter 
#' in a statistical model. May also be a list of lists of input matrices when a 
#' list of separate models is fit (e.g., with [flexsurvreg_list()]).
#' @param ... Arguments to pass to [id_attributes()].
#' @details Each row of each matrix `X` is an input vector, \eqn{x_{hik}}, where \eqn{h} denotes
#' a health-related index, \eqn{i} indexes a patient, and \eqn{k} is a treatment strategy. 
#' A health-related index is either a health state
#' (e.g., `state_id` or a transition between health states (e.g., `transition_id`).
#' In some cases, the health-related index \eqn{h} can be suppressed and separate models
#' can be fit for each health index. This is, for instance, the case in a [partitioned survival 
#' model][Psm] where separate models are fit for each survival endpoint. 
#' 
#' The rows of the matrices in `X` must be sorted in a manner consistent with the ID variables as
#' described in [id_attributes()].
#' @examples 
#' strategies <- data.frame(strategy_id = c(1, 2))
#' patients <- data.frame(patient_id = seq(1, 3), 
#'                           age = c(45, 47, 60),
#'                           female = c(1, 0, 0),
#'                           group = factor(c("Good", "Medium", "Poor")))
#' hesim_dat <- hesim_data(strategies = strategies,
#'                         patients = patients)
#' 
#' dat <- expand(hesim_dat, by = c("strategies", "patients"))
#' input_mats <- input_mats(X = list(mu = model.matrix(~ age, dat)),
#'                          strategy_id = dat$strategy_id,
#'                          n_strategies = length(unique(dat$strategy_id)),
#'                          patient_id = dat$patient_id,
#'                         n_patients = length(unique(dat$patient_id)))
#' print(input_mats)
#' @seealso [create_input_mats()]
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
check.input_mats <- function(object){
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

# Create input data from a fitted model ----------------------------------------
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

#' Check data argument for `create_input_mats`
#' 
#' Check that data argument for `create_input_mats` exists and that it is
#' of the correct type. 
#' @param data An object of class "expanded_hesim_data" returned by the function
#'  [expand.hesim_data]. 
#' @return If all tests passed, returns nothing; otherwise, throws an exception.
check_edata <- function(data){
  if (!inherits(data, "expanded_hesim_data")){
    stop("'data' must be of class 'expanded_hesim_data'.")
  } 
  if (!inherits(data, "data.table") & !inherits(data, "data.frame")){
    stop("'data' must inherit from either 'data.table' or 'data.frame'.")
  }   
  if (!inherits(data, "data.table")){
    setattr(data, "class", c("expanded_hesim_data", "data.table", "data.frame"))
  } 
}

extract_X <- function(coef_mat, data){
  varnames <- colnames(coef_mat)
  if (is.null(varnames)){
    stop("Variable names for coefficients cannot be NULL.")
  }
  if(!all(varnames %in% colnames(data))){
    stop("Not all variables in 'object' are contained in 'data'.",
         call. = FALSE)
  }  
  X <- as.matrix(data[, varnames, with = FALSE])
  return(X)
}

#' Create input matrices
#' 
#' \code{create_input_mats} is a generic function for creating an object of class
#' \code{\link{input_mats}}. Model matrices are typically constructed based on the 
#' variables specified in the model \code{object} and the data specified in \code{data}, 
#' although there are some cases in which \code{\link{input_mats}} can be created
#' from \code{object} alone.
#' @param object An object of the appropriate class. 
#' @param input_data An object of class "expanded_hesim_data" returned by the function
#'  \code{\link{expand.hesim_data}}. Used to look for the input variables needed to create an input matrix
#'  for use in a statistical models and the id variables for indexing rows in the input matrix. 
#' @param ... Further arguments passed to \code{\link{model.matrix}}.
#' @return An object of class \code{\link{input_mats}}.
#' @seealso \code{\link{input_mats}}.
#' @keywords internal
#' @examples 
#' library("flexsurv")
#' 
#' dt_strategies <- data.frame(strategy_id = c(1, 2))
#' dt_patients <- data.frame(patient_id = seq(1, 3), 
#'                           age = c(45, 47, 60),
#'                           female = c(1, 0, 0),
#'                           group = factor(c("Good", "Medium", "Poor")))
#' dt_states <- data.frame(state_id =  seq(1, 3),
#'                         state_name = factor(paste0("state", seq(1, 3))))
#' hesim_dat <- hesim_data(strategies = dt_strategies,
#'                         patients = dt_patients,
#'                         states = dt_states)
#'
#' # Class "lm"
#' expanded_dat <- expand(hesim_dat, by = c("strategies", "patients", "states"))
#' fit_lm <- stats::lm(costs ~ female + state_name, psm4_exdata$costs$medical)
#' input_mats <- create_input_mats(fit_lm, expanded_dat)
#' class(input_mats)
#'
#' # Class "flexsurvreg"
#' expanded_dat <- expand(hesim_dat, by = c("strategies", "patients"))
#' fit_wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                                  data = ovarian, dist = "weibull")
#' input_mats <- create_input_mats(fit_wei, expanded_dat)
#' class(input_mats)
#' @export
#' @rdname create_input_mats
create_input_mats <- function (object, ...) {
  if (missing(object)){
    stop("'object' is missing with no default.")
  }
  UseMethod("create_input_mats", object)
}

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
  check_edata(input_data)
  X_list <- formula_list_rec(object, input_data, ...)
  args <- c(list(X = X_list),
            get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

get_terms <- function(object){
  tt <- stats::terms(object)
  return(stats::delete.response(tt))
}

#' @export 
#' @rdname create_input_mats
create_input_mats.lm <- function(object, input_data, ...){
  check_edata(input_data)
  terms <- get_terms(object)
  X <- stats::model.matrix(terms, data = input_data, ...)
  args <- c(list(X = list(mu = X)),
           get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

#' @export 
#' @rdname create_input_mats
create_input_mats.lm_list <- function(object, input_data, ...){
  check_edata(input_data)
  X_list <- vector(mode = "list", length = length(object))
  names(X_list) <- names(object)
  for (i in 1:length(X_list)){
    terms <- get_terms(object[[i]])
    X_list[[i]] <- list(mu = stats::model.matrix(terms, data = input_data, ...))
  }
  args <- c(list(X = X_list),
           get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

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
  check_edata(input_data)
  X_list <- create_input_mats_flexsurvreg_X(object, input_data, ...)
  args <- c(list(X = X_list),
           get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

#' @export
#' @rdname create_input_mats
create_input_mats.flexsurvreg_list <- function(object, input_data,...){
  check_edata(input_data)
  X_list_2d <- vector(mode = "list", length = length(object))
  names(X_list_2d) <- names(object)
  for (i in 1:length(object)){
    X_list_2d[[i]] <- create_input_mats_flexsurvreg_X(object[[i]], input_data, ...)
  }
  args <- c(list(X = X_list_2d),
           get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

#' @export
#' @rdname create_input_mats
create_input_mats.partsurvfit <- function(object, input_data, ...){
  check_edata(input_data)
  return(create_input_mats.flexsurvreg_list(object$models, input_data, ...))
}

#' @export 
#' @rdname create_input_mats
create_input_mats.params_lm <- function(object, input_data, ...){
  check_edata(input_data)
  X <- extract_X(object$coefs, input_data)
  args <- c(list(X = list(mu = X)),
            get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

create_input_mats.params_surv_X <- function(object, input_data){
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
  check_edata(input_data)
  X_list <- create_input_mats.params_surv_X(object, input_data)
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
    X_list_2d[[i]] <- create_input_mats.params_surv_X(object[[i]], input_data)
  }
  args <- c(list(X = X_list_2d),
            get_input_mats_id_vars(input_data))
  return(do.call("new_input_mats", args))
}

create_input_mats_multinom_X <- function(object, input_data, ...){
  check_edata(input_data)
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
  args <- c(list(X = X ),
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
