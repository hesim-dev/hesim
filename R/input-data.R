# hesim data -------------------------------------------------------------------

#' Data table of treatment lines
#' 
#' Convert a list of treatment lines for multiple treatment strategies to a \code{\link{data.table}}
#' @param strategy_list A list where each element is a treatment strategy 
#' consisting of a vector of treatments. 
#' @param strategy_ids A numeric vector denoting the numeric id of each strategy
#' in \code{strategy_list}.
#' @return Returns a \code{\link{data.table}} in tidy format with three columns
#' \describe{
#' \item{strategy_id}{Treatment strategy ids.}
#' \item{line}{Line of therapy.}
#' \item{treatment_id}{Treatment id for treatment used at a given line of therapy within a treatment strategy.}
#' }
#' @examples 
#' strategies <- list(c(1, 2, 3),
#'                   c(1, 2))
#' lines_dt(strategies)
#' @export
lines_dt <- function(strategy_list, strategy_ids = NULL){
  treatments <- unlist(strategy_list, use.names = FALSE)
  if(!is.numeric(treatments)){
    stop("Elements in 'strategy_list' should be integers.")
  }
  n.treatments <- unlist(lapply(strategy_list, length), 
                                          use.names = FALSE)
  if (is.null(strategy_ids)){
    strategy_ids <- seq_len(length(strategy_list))
  }
  strategies <- rep(strategy_ids, times = n.treatments)
  lines <- unlist(lapply(n.treatments, seq_len), use.names = FALSE)
  return(data.table(strategy_id = strategies,
                    line = lines,
                    treatment_id = treatments))
}

#' Data for health-economic simulation modeling
#' 
#' A list of tables required for health-economic simulation modeling.
#' Each table must either be a \code{\link{data.frame}} or \code{\link{data.table}}. All id variables within 
#' each table must be numeric vectors of integers. 
#' @param strategies A table of treatment strategies. 
#' Must contain the column \code{strategy_id} denoting a unique strategy. Other columns are variables
#'  describing the characteristics of a treatment strategy. 
#' @param patients A table of patient observations. 
#' Must contain the column \code{patient_id} denoting a unique patient. The 
#' number of rows should be equal to the number of patients in the model.
#' Other columns are variables describing the characteristics of a patient.
#' @param lines A table of treatment lines used for each treatment strategy. Must contain the columns
#' \code{strategy_id}, denoting a treatment strategy, and \code{line}, denoting a treatment line. Other 
#' columns are variables describing the characteristics of a treatmnet line for a given treatment
#' strategy. A column denoting the treatment used for a given strategy and line would often
#' be specified.
#' @param states A table of health states. Must contain the column
#' \code{state_id}, which denotes a unique health state. The number of rows should
#' be equal to the number of health states in the model. Other columns can describe the
#' characteristics of a health state.
#' @param transitions A table of health state transitions. Must contain the column
#' \code{transition_id}, which denotes a unique transition; \code{from}, which denotes
#' the starting health state; and \code{to} which denotes the state that will be
#' transitioned to.
#' @return Returns an object of class "hesim_data", which is a list of data tables for
#' health economic simulation modeling.
#' @examples 
#' dt.strategies <- data.frame(strategy_id = c(1, 2))
#' dt.patients <- data.frame(patient_id = seq(1, 3), age = c(65, 50, 75),
#'                           gender = c("Female", "Female", "Male"))
#' dt.lines <- lines_dt(list(c(1, 2, 5), c(1, 2)))
#' dt.states <- data.frame(state_id =  seq(1, 3),
#'                         state_var = c(2, 1, 9))
#' hesim.dat <- hesim_data(strategies = dt.strategies,
#'                           patients = dt.patients,
#'                           states = dt.states,
#'                           lines = dt.lines)
#' @export
hesim_data <- function(strategies, patients, lines = NULL, states = NULL,
                          transitions = NULL){
  object <- new_hesim_data(strategies = strategies, patients = patients,
                           lines = lines, states = states, transitions = transitions)
  return(check(object))
}

new_hesim_data <- function(strategies, patients, lines = NULL, states = NULL,
                          transitions = NULL){
  data <- list()
  data$strategies <- strategies 
  data$patients <- patients
  data$lines <- lines
  data$states <- states
  data$transitions <- transitions
  class(data) <- c("hesim_data")
  return(data)
}

check.hesim_data <- function(x){
  # Strategies
  if (is.null(x$strategies)){
    stop("'strategies' cannot be NULL'.",
         call. = FALSE)
  }
  check_hesim_data_type(x$strategies, "strategies")
  if (!"strategy_id" %in% colnames(x$strategies)){
    stop("'strategy_id' must be a column of 'strategies'.",
         call. = FALSE)
  }
  
  # Patients
  if (!is.null(x$patients)){
      check_hesim_data_type(x$patients, "patients")
      if (!"patient_id" %in% colnames(x$patients)){
        stop("'patient_id' must be a column of 'patients'.", 
             call. = FALSE)
      }
  }
  
  # Lines
  if (!is.null(x$lines)){
      check_hesim_data_type(x$lines, "lines")
      if (!"strategy_id" %in% colnames(x$lines)){
        stop("'strategy_id' must be a column of 'lines'.", 
             call. = FALSE)
      }
      if (!"line" %in% colnames(x$lines)){
        stop("'line' must be a column of 'lines'.", 
             call. = FALSE)
      }
  }
  
  # States
  if (!is.null(x$states)){
      check_hesim_data_type(x$states, "states")
      if (!"state_id" %in% colnames(x$states)){
        stop("'state_id' must be a column of 'states'.", 
             call. = FALSE)
      }
  }
  
  # Transitions
  if (!is.null(x$transitions)){
      check_hesim_data_type(x$transitions, "transitions")
      if (!"transition_id" %in% colnames(x$transitions)){
        stop("'transition_id' must be a column of 'transitions'.", 
             call. = FALSE)
      }
  }

  return(x)
}

check_hesim_data_type <- function(tbl, tbl_name){
  if(!is.data.frame(tbl) & !is.data.table(tbl)){
    msg <- paste0("'", tbl_name, "'", " must be a data.frame or data.table.")
    stop(msg, call. = FALSE)
  }
}

#' Expand hesim_data
#' 
#' Expand the data tables from an object of class \code{\link{hesim_data}} into a single \code{\link{data.table}} 
#' of class "expanded_hesim_data". See the description of the return value for details on how this is done.
#' @param data An object of class \code{\link{hesim_data}}.
#' @param by A character vector of the names of the data tables in \code{\link{hesim_data}} to expand by.
#' @return This function is similar to \code{\link{expand.grid}}, but works for data frames or data tables. 
#' Specifically, it creates a \code{data.table} from all combinations of the supplied tables in \code{data}. 
#' The supplied tables are determined using the \code{by} argument. The resulting dataset is sorted by 
#' prioritizing id variables as follows: (i) \code{strategy_id}, (ii) \code{line}, (iii) \code{patient_id},
#' (iv) the health-related id variable (either \code{state_id} or \code{transition_id}).
#' @examples 
#' dt.strategies <- data.frame(strategy_id = c(1, 2))
#' dt.patients <- data.frame(patient_id = seq(1, 3), age = c(65, 50, 75),
#'                           gender = c("Female", "Female", "Male"))
#' dt.lines <- lines_dt(list(c(1, 2, 5), c(1, 2)))
#' dt.states <- data.frame(state_id =  seq(1, 3),
#'                         state_var = c(2, 1, 9))
#' data.tables <- hesim_data(strategies = dt.strategies,
#'                           patients = dt.patients,
#'                           states = dt.states,
#'                           lines = dt.lines)
#' expand_hesim_data(data.tables, by = c("strategies", "patients"))
#' @export
expand_hesim_data <- function(data, by = c("strategies", "patients")){
  if ("transitions" %in% by & "states" %in% by){
    stop("Cannot expand by both 'transitions' and 'states'.", call. = FALSE)
  }
  if (!all(by %in% names(hesim_data_sorting_map()))){
    stop("One of the elements specified in 'by' is not a table in 'hesim_data'.",
         call. = FALSE)
  }
  sorted.by <- hesim_data_sorted_by(by)
  tbl.list <- data[sorted.by]
  for (i in 1:length(tbl.list)){
    if (is.null(tbl.list[[i]])){
      stop("Cannot merge a NULL data table.")
    }
  }
  tbl.list <- lapply(tbl.list, function(x) as.data.frame(x))
  dat <- if ("lines" %in% sorted.by){
    tbl.list1 <- tbl.list[which(names(tbl.list) != "lines")]
    dat <- Reduce(function(...) merge(..., by = NULL), tbl.list1)
    dat <- data.table(merge(dat, tbl.list[[which(names(tbl.list) == "lines")]],
                     by = "strategy_id", sort = FALSE))
  } else{
    dat <- data.table(Reduce(function(...) merge(..., by = NULL, sort = FALSE), tbl.list))
  }
  sort_hesim_data(dat, sorted.by)
  id.cols <- unlist(hesim_data_sorting_map()[sorted.by])
  nonid.cols <- colnames(dat)[!colnames(dat) %in% id.cols]
  dat <- dat[, c(id.cols, nonid.cols), with = FALSE]
  res <- list(data = dat,
              id_vars = unname(id.cols))
  class(res) <- "expanded_hesim_data"
  return(res)
}

hesim_data_sorting_map <- function(){
  list(strategies = "strategy_id",
       lines = "line",
       patients = "patient_id",
       transitions = "transition_id",
       states = "state_id")
}

hesim_data_sorted_by <- function(by){
  sorted.map <- hesim_data_sorting_map()
  sorted.by <- names(sorted.map)[names(sorted.map) %in% by]
  return(sorted.by)
}

sort_hesim_data <- function(data, sorted_by){
  setorderv(data, unlist(hesim_data_sorting_map()[sorted_by]))
}


# input_data class -------------------------------------------------------------
#' Input data for a statistical model
#' 
#' Create an object of class "input_data", which contains data for use
#' as an input to a statistical model. Data consists of (i) input matrices, \code{X},
#' (ii) id variables indexing the rows of each matrix in \code{X}, and (iii) the dimensions of
#' the \code{X} matrices. More details are provided under "Details" below. Note that an "input_data" 
#' object can often be created more easily using the generic function 
#' \code{\link{form_input_data}}. 
#' 
#' @param X A list of input matrices for predicting the values of each parameter in a statistical model. May also be
#' a list of lists of input matrices when a list of separate models is fit (e.g., with \link{flexsurvreg_list}).
#' @param strategy_id A numeric vector of integers denoting the treatment strategy represented by each row
#' in \code{X}.
#' @param n_strategies A scalar denoting the number of unique treatment strategies.
#' @param patient_id A numeric vector of integers denoting the patient represented by each row
#' in \code{X}.
#' @param n_patients A scalar denoting the number of unique patients.
#' @param line A numeric vector of integers denoting the treatment line represented by each row
#' in \code{X}. Not supported by currently available models.
#' @param n_lines A \code{\link{data.table}} denoting the number of treatment lines associated
#' with each treatment strategy. Should contain a column, "strategy_id", and a column,
#' "N". Not supported by currently available models.
#' @param state_id A numeric vector of integers denoting the health state represented by each row
#' in \code{X}.
#' @param n_states A scalar denoting the number of unique health states.
#' @param transition_id A numeric vector denoting the 
#' health state transition represented by each row in \code{X}. This must only be specified when
#' estimating the health state transitions with a joint likelihood function. If independent
#' models are fit for each transition, then separate \code{X} matrices must be specified
#' for each transition. Note that this is not currently supported but
#' will be supported once \code{hesim} provides support for state transition modeling.
#' @param n_transitions A scalar denoting the number of unique transitions. 
#' Not supported by currently available models. 
#' @param time_fun A pointer to a C++ functor that can be used to update \code{X} as a function
#' of time in a simulation model. Not currently supported.
#' 
#' @details Each row of each matrix \code{X} is an input vector, \eqn{x_{hijk}}, where \eqn{h} denotes
#' a health-related index, \eqn{i} indexes a patient, \eqn{j} indexes a treatment line,
#' and \eqn{k} is a treatment strategy. A health-related index is either a health state
#' (e.g., \code{state_id}) or a transition between health states (e.g., \code{transition_id}).
#' In some cases, the health-related index \eqn{h} can be suppressed and separate models
#' can be fit for each health index. This is, for instance, the case in a partitioned survival 
#' model where separate models are fit for each survival endpoint. Likewise, models 
#' can be fit without multiple treatment lines as would, again, be the case in a 
#' partitioned survival analysis where sequential treatment would be incorporated by 
#' adding additional health states rather than by using the index \eqn{j}.
#' 
#' The rows of the matrices in \code{X} must be sorted in a manner consistent with the id variables.
#' The sorting order should be the same as specified in \code{\link{expand_hesim_data}}; that is,
#' the rows of \code{X} must be sorted by: (i) \code{strategy_id}, (ii) \code{line}, 
#' (iii) \code{patient_id}, and (iv) the health-related id variable (either \code{state_id} or
#'  \code{transition_id}).
#' @examples 
#' dt.strategies <- data.frame(strategy_id = c(1, 2))
#' dt.patients <- data.frame(patient_id = seq(1, 3), 
#'                           age = c(45, 47, 60),
#'                           female = c(1, 0, 0),
#'                           group = factor(c("Good", "Medium", "Poor")))
#' hesim.dat <- hesim_data(strategies = dt.strategies,
#'                              patients = dt.patients)
#' 
#' dat <- expand_hesim_data(hesim.dat, by = c("strategies", "patients"))$data
#' input.dat <- input_data(X = list(mu = model.matrix(~ age, dat)),
#'                        strategy_id = dat$strategy_id,
#'                        n_strategies = length(unique(dat$strategy_id)),
#'                        patient_id = dat$patient_id,
#'                        n_patients = length(unique(dat$patient_id)))
#' print(input.dat)
#' @export
input_data <- function(X, strategy_id, n_strategies,
                       patient_id, n_patients,
                       line = NULL, n_lines = NULL,
                       state_id = NULL, n_states = NULL,
                       transition_id = NULL, n_transitions = NULL,
                       time_fun = NULL){
  object <- new_input_data(X, strategy_id, n_strategies,
                                    patient_id, n_patients,
                                    line, n_lines,
                                    state_id, n_states,
                                    transition_id, n_transitions,
                                    time_fun)
  check(object)
  return(object)
}

new_input_data <- function(X, strategy_id, n_strategies,
                       patient_id, n_patients,
                       line = NULL, n_lines = NULL,
                       state_id = NULL, n_states = NULL,
                       transition_id = NULL, n_transitions = NULL,
                       time_fun = NULL){
  stopifnot(is.matrix(X) | is.list(X))
  stopifnot(is.numeric(strategy_id))
  stopifnot(is.numeric(n_strategies))
  stopifnot(is.numeric(patient_id))
  stopifnot(is.numeric(n_patients))
  stopifnot(is.numeric(line) | is.null(line))
  stopifnot(is.data.table(n_lines) | is.null(n_lines))
  stopifnot(is.numeric(state_id) | is.null(state_id))
  stopifnot(is.numeric(n_states) | is.null(n_states))
  stopifnot(is.numeric(transition_id) | is.null(n_transitions))
  stopifnot(is.numeric(n_transitions) | is.null(n_transitions))
  if(!is.null(time_fun)){
   if(!inherits(time_fun, "externalptr")){
     stop("If not NULL, 'time_fun' must be of class 'externalptr'.",
          call. = FALSE)
   } 
  }
  
  object <- list(X = X,
                 strategy_id = strategy_id, n_strategies = n_strategies,
                 patient_id = patient_id, n_patients = n_patients,
                 line = line, n_lines = n_lines,
                 state_id = state_id, n_states = n_states,
                 transition_id = transition_id, n_transitions = n_transitions,
                time_fun = time_fun)
  object[sapply(object, is.null)] <- NULL
  class(object) <- "input_data"
  return(object)
}

#' @rdname check
check.input_data <- function(object){
  if (!is.list(object$X)){
    stop("'X' must be a list or a list of lists.", call. = FALSE)
  }
  X <- flatten_lists(object$X)
  if(!all(sapply(X, function(y) inherits(y, "matrix")))){
    stop("'X' must be a list or list of lists of matrices.", call. = FALSE)
  }
  X.nrows <- sapply(X, nrow)
  if (is.list(object$X)){
    if (!all(X.nrows[1] == X.nrows)){
      stop("The number of rows in each matrix in 'X' must be the same.",
           call. = FALSE)
    }
  }
  
  # Check that strategies in n_lines are correct
  if (!is.null(object$n_lines)){
    if(!all(sort(object$n_lines$strategy_id) == sort(unique(object$strategy_id)))){
      msg <- "'strategy_id' in 'n_lines' is not consistent with 'strategy_id'."
      stop(msg, call. = FALSE)
    }
  }
  
  id.vars <- c("strategy_id", "patient_id", "line", "state_id", "transition_id")
  id.vars.n <- c("n_strategies", "n_patients", "n_lines", "n_states", "n_transitions")
  for (i in 1:length(id.vars)){
    if (!is.null(object[[id.vars[i]]])){
      
      ## Check number of rows
      if(length(object[[id.vars[i]]]) != X.nrows[1]){
        msg <- paste0("The length of '", id.vars[i], "' does not equal the number of rows in the ",
                      "'X' matrices.")
        stop(msg, call. = FALSE)
      }
      
      # Check that n_strategies, n_patients, ..., is correct
        if (id.vars[i] != "line"){
          if(length(unique(object[[id.vars[i]]])) != object[id.vars.n[i]]){
            msg <- paste0("The number of unique observations in '", id.vars[i], 
                      "' does not equal '", id.vars.n[i], "'.")
            stop(msg, call. = FALSE)
          } 
          } else{
            line <- object[[id.vars[i]]]
            expected.n <- data.table("strategy_id" = object[["strategy_id"]],
                                     "line" = line)
            expected.n <- expected.n[, list(N = length(unique(line))), by = "strategy_id"]
            if(!all(expected.n$N == object[[id.vars.n[i]]]$N)){
              msg <- paste0("The number of unique observations in 'line'", 
                        "' within each 'strategy_id' does not equal 'n_lines$N'")
              stop(msg, call. = FALSE)
            }  
        } # end if else for line id variable
    } # end loop of id.vars
  }
  
  # Check if id variables are sorted properly 
  indices.df <- data.table(do.call("cbind", object[id.vars]))
  sorted.seq <- seq_len(nrow(indices.df))
  indices.df[, "row_num" := sorted.seq]
  by <- id.vars[sapply(object[id.vars], function(x) !is.null(x))]
  sort_hesim_data(indices.df, sorted_by = hesim_data_sorted_by(by))
  if(!all(indices.df$row_num == sorted.seq)){
    msg <- paste0("The id variables are not sorted correctly. The sort priority of the ",
                  "id variables must be as follows: 'strategy_id', 'line', 'patient_id' and ",
                   "the health-related id variable ('state_id' or 'transition_id').")
    stop(msg, call. = FALSE)
  }

  # Check if the number of unique observations is correct within groups
  indices.df[, "row_num" := NULL]
  for (i in 2:ncol(indices.df)){
    dt.by <- colnames(indices.df)[i - 1]
    col <- colnames(indices.df)[i]
    if (id.vars.n[i] == "n_lines"){
      len <- indices.df[, list(len = length(unique(get(col)))), 
                        by = c("strategy_id", dt.by)]$len
    } else{
      len <- indices.df[, list(len = length(unique(get(col)))), by = dt.by]$len 
    }
    if (id.vars.n[i] == "n_lines"){
      user.n <- object[[id.vars.n[i]]]$N
    } else{
      user.n <- object[[id.vars.n[i]]]
    }
    if (!all(unique(len) == user.n)){
      if (id.vars.n[i] == "n_lines"){
        msg <- paste0("The number of unique '", col, "' observations within each value",
                    " of 'strategy_id' and '", dt.by, " ' must be consistent with '", id.vars.n[i], "'.")
        stop(msg, call. = FALSE)
      } else{
        msg <- paste0("The number of unique '", col, "' observations within each value",
                    " of '", dt.by, " ' must equal '", id.vars.n[i], "'.")
        stop(msg, call. = FALSE)
      }
    }
  }
}

# Form input data from a fitted model ------------------------------------------
size_id_map <- function(){
  c(strategy_id = "n_strategies", 
    patient_id = "n_patients",
    line = "n_lines",
    state_id = "n_states",
    transition_id = "n_transitions")
}

get_input_data_id_vars <- function(data){
  map <- size_id_map()
  res <- list() 
  id.vars <- data$id_vars
  for (i in 1:length(id.vars)){
    res[[id.vars[i]]] <- data$data[[id.vars[i]]]
    if (id.vars[i] != "line"){
       res[[map[id.vars[i]]]] <- length(unique(data$data[[id.vars[i]]]))
    } else{
      n.lines <- data$data[, .N, by = c("strategy_id", "line")][, .N, by = "strategy_id"]
      res[[map[id.vars[i]]]] <- n.lines
    }
  }
  return(res)
}

#' Form input data
#' 
#' \code{form_input_data} is a generic function for forming data that can be
#' used as a input to a statistical model. Model matrices are formed based on the 
#' variables specified in the model \code{object} and the data specified in \code{data}.
#' @param object An object of the appropriate class. Currently supports
#' \code{\link{formula_list}}, \code{\link{lm}}, \code{\link{flexsurvreg}}, 
#'  \code{\link{flexsurvreg_list}}, and \code{\link{partsurvfit}}.
#' @param data An object of class "expanded_hesim_data" returned by the function
#'  \code{\link{expand_hesim_data}}. Used to look for the input variables needed to create an input matrix
#'  for use in a statistical models and the id variables for indexing rows in the input matrix. 
#' @param ... Further arguments passed to \code{\link{model.matrix}}.
#' @return An object of class \code{\link{input_data}}.
#' @examples 
#' library("flexsurv")
#' 
#' dt.strategies <- data.frame(strategy_id = c(1, 2))
#' dt.patients <- data.frame(patient_id = seq(1, 3), 
#'                           age = c(45, 47, 60),
#'                           female = c(1, 0, 0),
#'                           group = factor(c("Good", "Medium", "Poor")))
#' dt.states <- data.frame(state_id =  seq(1, 3),
#'                         state_name = factor(paste0("state", seq(1, 3))))
#' hesim.dat <- hesim_data(strategies = dt.strategies,
#'                         patients = dt.patients,
#'                         states = dt.states)
#'
#' # Class "lm"
#' expanded.dat <- expand_hesim_data(hesim.dat, by = c("strategies", "patients", "states"))
#' fit.lm <- stats::lm(costs ~ female + state_name, part_surv4_simdata$costs$medical)
#' input.dat <- form_input_data(fit.lm, expanded.dat)
#' class(input.dat)
#'
#' # Class "flexsurvreg"
#' expanded.dat <- expand_hesim_data(hesim.dat, by = c("strategies", "patients"))
#' fit.wei <- flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1, 
#'                                  data = ovarian, dist = "weibull")
#' input.dat <- form_input_data(fit.wei, expanded.dat)
#' class(input.dat)
#' @export
form_input_data <- function (object, data, ...) {
  if (missing(object)){
    stop("'object' is missing with no default.")
  }
  if (missing(data)){
    stop("'data' is missing with no default.")
  }
  if (!inherits(data, "expanded_hesim_data")){
    stop("'data' must be of class 'expanded_hesim_data'.")
  }
  UseMethod("form_input_data", object)
}

formula_list_rec <- function(object, data, ...){
  x <- vector(mode = "list", length = length(object))
  names(x) <- names(object)
  for (i in 1:length(x)){
    if (inherits(object[[i]], "formula")){
      x[[i]] <- stats::model.matrix(object[[i]], data = data$data, ...)
    } else{
      x[[i]] <- formula_list_rec(object[[i]], data = data, ...)
    }
  }
  return(x)
}

#' @export
form_input_data.formula_list <- function(object, data, ...){
  X.list <- formula_list_rec(object, data, ...)
  args <- c(list(X = X.list),
           get_input_data_id_vars(data))
  return(do.call("new_input_data", args))
}

get_terms <- function(object){
  tt <- stats::terms(object)
  return(stats::delete.response(tt))
}

#' @export 
form_input_data.lm <- function(object, data, ...){
  terms <- get_terms(object)
  X <- stats::model.matrix(terms, data = data$data, ...)
  args <- c(list(X = list(mu = X)),
           get_input_data_id_vars(data))
  return(do.call("new_input_data", args))
}

#' @export 
form_input_data.lm_list <- function(object, data, ...){
  X.list <- vector(mode = "list", length = length(object))
  names(X.list) <- names(object)
  for (i in 1:length(X.list)){
    terms <- get_terms(object[[i]])
    X.list[[i]] <- list(mu = stats::model.matrix(terms, data = data$data, ...))
  }
  args <- c(list(X = X.list),
           get_input_data_id_vars(data))
  return(do.call("new_input_data", args))
}

form_input_data_flexsurvreg_X <- function(object, data, ...){
  pars <- object$dlist$pars
  X.list <- vector(mode = "list", length = length(pars))
  names(X.list) <- pars
  for (i in 1:length(pars)){
    form <- object$all.formulae[[pars[i]]]
    if (is.null(form)){
      form <- stats::formula(~1)
    } else{
      form <- stats::delete.response(stats::terms(form))
    }
    X.list[[i]] <- stats::model.matrix(form, data = data$data, ...)
  }
  return(X.list)
}

#' @export
form_input_data.flexsurvreg <- function(object, data,...){
  X.list <- form_input_data_flexsurvreg_X(object, data, ...)
  args <- c(list(X = X.list),
           get_input_data_id_vars(data))
  return(do.call("new_input_data", args))
}

#' @export
form_input_data.flexsurvreg_list <- function(object, data,...){
  X.list.2d <- vector(mode = "list", length = length(object))
  names(X.list.2d) <- names(object)
  for (i in 1:length(object)){
    X.list.2d[[i]] <- form_input_data_flexsurvreg_X(object[[i]], data, ...)
  }
  args <- c(list(X = X.list.2d),
           get_input_data_id_vars(data))
  return(do.call("new_input_data", args))
}

#' @export
form_input_data.partsurvfit <- function(object, data, ...){
  return(form_input_data.flexsurvreg_list(object$models, data, ...))
}

#' @export
form_input_data.joined_flexsurvreg_list <- function(object, data,...){
  models <- object$models
  X.list.3d <- vector(mode = "list", length = length(models))
  names(X.list.3d) <- names(models)
  for (i in 1:length(models)){
    X.list.3d[[i]] <- vector(mode = "list", length = length(models[[i]]))
    names(X.list.3d[[i]]) <- names(models[[i]])
    for (j in 1:length(models[[i]])){
      X.list.3d[[i]][[j]] <- form_input_data_flexsurvreg_X(models[[i]][[j]], data, ...)
    }
  }
  args <- c(list(X = X.list.3d),
           get_input_data_id_vars(data))
  return(do.call("new_input_data", args))
}
